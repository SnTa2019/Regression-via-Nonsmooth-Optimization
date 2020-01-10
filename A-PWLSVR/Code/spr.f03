!*************************************************************************
!*                                                                       *
!*     A-PWLSVR - Adaptive Piecewise Linear Support Vector Regression    *
!*                                                                       *
!*     by Napsu Karmitsa 2019, Adil Bagirov and Kaisa Joki               *
!*     (last modified 08.01.2020).                                       *
!*                                                                       *
!*     The software is free for academic teaching and research           *
!*     purposes but we ask you to refer the appropriate references       *
!*     given below, if you use it.                                       *
!*                                                                       *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*
!*     spr.f03               - Mainprogram (this file).
!*     constants.f03         - Parameters and constants.
!*     initspr.f03           - Initialization of parameters for A-PWLSVR and DBDC.
!*     functions.f03         - Computation of DC components f_1 and f_2 and their subgradients for the PWLSVR problem.
!*     bundle1.f03           - Bundle of DC component f_1.
!*     bundle2.f03           - Bundle of DC component f_2.
!*     dbdc.f03              - DBDC method.
!*     plqdf1.f              - Quadratic solver by L. Luksan.
!*
!*     Makefile              - makefile.
!*
!*
!*     Calling sequence:
!*
!*       ./spr infile nrecord ntrain nft maxmin [noutcom maxlin eps]
!*
!*     with the following arguments:
!*
!*     infile                - the data set;
!*     nrecord               - number of records in the data set;
!*     ntrain                - size of the training set, ntrain <= nrecord;
!*     nft                   - number of features in the data set;
!*     maxmin                - maximum number of min functions under maximum in the model;
!*     noutcom               - number of output column, optional, nft by default;
!*     maxlin                - number of linear functions under each minimum in the model, optional,
!*                             Default value = 3;
!*     eps                   - stopping tolerance, optional, 10^(-2) by default.
!*
!*     If you want to tune other parameters, modify initspr.f03. 
!*     As default the results are written in
!*        
!*     predictions.txt       - Evaluation criteria (e.g. RMSE, MAE,...) with training and predictions;
!*     reg_results.txt       - Regression results with hyperplanes.
!*     pwl_points.txt        - Points of original data and pwl-model obtained
!*                             for drawing purposis. Works only with nft = 2 and maxlin = 3 (otherwise empty). 
!*
!*
!*     References:
!*       A. Bagirov, S. Taheri, N. Karmitsa, K. Joki, M.M. Mäkelä, "Adaptive piecewise linear support vector regression", submitted 2020.
!*
!*
!*     Acknowledgements:
!*
!*       The work was financially supported by the Academy of Finland (Project No. 289500 and 319274).
!*
!*
!*****************************************************************************************
!*
!*     * PROGRAM spr *
!*
!*     Main program for A-PWLSVR piecewise linear support vector regression.
!*
!*****************************************************************************************


PROGRAM spr

    USE constants, ONLY : dp   ! Precision for reals.         
    
    USE initspr, ONLY : &      ! Initialization of parameters.
        infile, &              ! Name of the dataset file.
        outfile1, &            ! Regression results with hyperplanes.
        outfile2, &            ! Regression results with predictions.
        outfile3, &            ! Data points and pwl-model (only usable for two dim data).
        maxv, &                ! Maximum value of data points.
        minv, &                ! Minimum value of data points.
        nrecord, &             ! Number of instances in data, from user.
        ntrain, &              ! Size of the training set, ntrain <= nrecord, from user.
        nft, &                 ! Number of features in data, from user.
        noutcom, &             ! Output column.
        iscale, &              ! Scaling is
                               !   0 - not used,
                               !   1 - used.
        inc, &                 ! Incremental algorithm is
                               !   0 - not used,
                               !   1 - used.
        ipredi, &              ! Prediction is computed
                               !   0 - only at the termination with complite model,
                               !   1 - at every step.                                
        maxmin, &              ! Maximum number of min functions under maximum in the model.
        maxmin2, &             ! maxmin2 = 1,...,maxmin.
        maxlin, &              ! Number of linear functions under each minimum in the model.
        mk, &                  ! Vector of the number of linear functions under each
                               ! minimum functions (for the future, now always the same number).
        maxdim, &              ! Maximum dimension of the problem.                  
        tol, &                 ! Stopping tolerance.
        eps, &                 ! Epsilon used in SVM.
        a, &                   ! Data matrix a(nrecord,nft), from input file.
        a4, &                  ! Auxiliary data matrix a4(ntrain,nft).
        x_0=>x0, &             ! The starting point. 
        x0reg, &               ! The SVR solution.
        x_solution, &          ! The solution.
        x_solution2, &         ! The solution.
        x_unscaled, &          ! The unscaled solution.
        f_solution, &          ! The objective function value at the solution 'x_solution'
        f_solution2, &         ! The objective function value at the solution 'x_solution2'
        pred, &                ! the values of predicted outputs
        dlin, &
        clin, & 
        init_par, &
        scaling, &
        rescaling
  
    USE initdbdc, ONLY : &     ! Initialization of parameters.
        problem1, &
        problem2, &  
        agg_used, &
        stepsize_used, &
        termination, &
        counter, &
        iprint, &
        mit, &
        mrounds, &
        mrounds_clarke
         
         
    USE bundle1                ! The BUNDLE of the DC component f_1
    USE bundle2                ! The BUNDLE of the DC component f_2
    USE dbdc_non               ! DBDC method


    IMPLICIT NONE

    !local
    REAL(KIND=dp) :: &
        h,h2,tmp,tmp2, &
        atmp, &
        rmse1, &
        mae1, &
        mse1, &
        ce1, &
        r1, &
        e0_test,e0_tr,e1_test,e1_tr, &
        f_1,f_2, &
        sx,sy,rcor,ytest, &
        large,piece,piece_best,piece_min, &
        start_time, finish_time   ! start and finish CPU time
    INTEGER :: &
        count, &               ! Number of commandline arguments
        ntest, &               ! Size of the test set.
        neval, &               ! Total number of linear function evaluations.
        i,i2,i3,j,j1,j2,j3, &
        ipiecemin,ipiecebest, &
        pind, &
        mtmp
        
    
    CHARACTER(len=50) :: arg

    ! Intrinsic Functions
    INTRINSIC ABS,MAX

    count = command_argument_count()
    IF (count < 1) THEN
        PRINT*
        PRINT*,'A-PWLSVR     - Adaptive Piecewise Linear Support Vector Regression Method'
        PRINT*
        PRINT*,'Calling sequence:'
        PRINT*
        PRINT*,'   ./spr infile nrecord ntrain nft maxmin [noutcom maxlin tol]'
        PRINT*
        PRINT*,'with the following arguments:'
        PRINT*
        PRINT*,'   infile    - the data set;'
        PRINT*,'   nrecord   - number of records in the data set;'
        PRINT*,'   ntrain    - number of records in the training set;'
        PRINT*,'   nft       - number of features in the data set;'
        PRINT*,'   maxmin    - maximum number of min functions under maximum in the model;'
        PRINT*,'   noutcom   - optional, number of output column, nft by default;'
        PRINT*,'   maxlin    - optional, number of linear functions under each minimum in the model,'
        PRINT*,'               2 by default;'
        PRINT*,'   tol       - optional, stopping tolerance, 10^(-2) by default.'
        PRINT*
        PRINT*,'If you want to tune other parameters, modify initspr.f03.' 
        PRINT*,'As default the results are written in'
        PRINT*
        PRINT*,'   predictions.txt  - Evaluation criteria (e.g. RMSE, MAE,...) with training and predictions;'
        PRINT*,'   reg_results.txt  - Regression results with hyperplanes.'
        PRINT*,'   pwl_points.txt   - Points of original data and pwl-model obtained for drawing purposis.'
        PRINT*,'                      Works only with nft = 2 and maxlin = 3 (otherwise empty).' 
        PRINT*
        STOP

    END IF

    IF (count < 5) STOP '  Too few commandline arguments. Type ./spr for help.'

    CALL get_command_argument(2, arg)
    READ (arg,*) nrecord
    CALL get_command_argument(3, arg)
    READ (arg,*) ntrain
    CALL get_command_argument(4, arg)
    READ (arg,*) nft
    CALL get_command_argument(5, arg)
    READ (arg,*) maxmin

    IF (count >= 6) THEN
        CALL get_command_argument(6, arg)
        READ (arg,*) noutcom
    ELSE
        noutcom=nft    
    END IF
    IF (count >= 7) THEN
        CALL get_command_argument(7, arg)
        READ (arg,*) maxlin
    END IF
    IF (count >= 8) THEN
        CALL get_command_argument(8, arg)
        READ (arg,*) tol
    END IF


    CALL get_command_argument(1, infile)
    
    PRINT*,' Input file: ',infile
    PRINT*,' Number of records = ',nrecord
    PRINT*,' Size of the training set = ',ntrain
    PRINT*,' Number of features = ',nft
    PRINT*,' Maximum number of min functions = ',maxmin
    PRINT*,' Maximum number of linear functions = ',maxlin

    CALL init_par() ! sets dimension for optimization problem etc.
    
    allocate(a(nrecord,nft),a4(ntrain,nft))
    allocate(x_0(maxdim),x0reg(maxdim),x_solution(maxdim),x_solution2(maxdim),x_unscaled(maxdim),clin(nft),dlin(nft)) ! maxdim= nft*q, where q=sum_k=1,...K (M_k)
    allocate(pred(MAX(nrecord-ntrain,ntrain)))
    allocate(mk(maxmin))

    neval=0
    large=1.0E+36_dp    
    mk=maxlin ! each min function has maxlin linear functions
          
    OPEN(43,file=outfile3)
    OPEN(42,file=outfile1)
    OPEN(41,file=outfile2)
    OPEN(78,file=infile,status='old',form='formatted')

    DO i=1,nrecord
        READ(78,*) (a(i,j),j=1,nft)
    END DO
        
    SELECT CASE (inc)

        CASE(0)
            WRITE(42,*) 'Piecewise linear regression without incremental approach.'
            WRITE(42,*) 'Data set: ',infile
            WRITE(42,100) maxmin,maxlin,eps
            IF(iscale == 1) WRITE(42,*) 'Scaling is used. '
            IF(iscale == 0) WRITE(42,*) 'No scaling. '
            
        
            WRITE(41,*) 'Piecewise linear regression without incremental approach.'
            WRITE(41,*) 'Data set: ',infile
            WRITE(41,100) maxmin,maxlin,eps
            IF(iscale == 1) WRITE(41,*) 'Scaling is used. '
            IF(iscale == 0) WRITE(41,*) 'No scaling. '
            
        
        CASE(1)
            WRITE(42,*) 'A-PWLSVR - Adaptive Piecewise linear support vector regression.'
            WRITE(42,*) 'Data set: ',infile
            WRITE(42,100) maxmin,maxlin,eps
            IF(iscale == 1) WRITE(42,*) 'Scaling is used. '
            IF(iscale == 0) WRITE(42,*) 'No scaling. '
            
        
            WRITE(41, *) 'A-PWLSVR - Adaptive Piecewise linear support vector regression.'
            WRITE(41, *) 'Data set: ',infile
            WRITE(41,100) maxmin,maxlin,eps
            IF(iscale == 1) WRITE(41,*) 'Scaling is used. '
            IF(iscale == 0) WRITE(41,*) 'No scaling. '
            
        CASE(2)
            WRITE(42,*) 'A-PWLSVR - Adaptive Piecewise linear support vector regression with double starting points.'
            WRITE(42,*) 'Data set: ',infile
            WRITE(42,100) maxmin,maxlin,eps
            IF(iscale == 1) WRITE(42,*) 'Scaling is used. '
            IF(iscale == 0) WRITE(42,*) 'No scaling. '


            WRITE(41, *) 'A-PWLSVR - Adaptive Piecewise linear support vector regression with double starting points.'
            WRITE(41, *) 'Data set: ',infile
            WRITE(41,100) maxmin,maxlin,eps
            IF(iscale == 1) WRITE(41,*) 'Scaling is used. '
            IF(iscale == 0) WRITE(41,*) 'No scaling. '
        CASE DEFAULT
            
    END SELECT
    
100 FORMAT(' Model with maxmin = ',i3,', maxlin = ', i3,', and epsilon = ',f5.3)
    

        
    IF (noutcom < nft) THEN ! swap features

        DO i = 1,nrecord
            atmp = a(i,noutcom)
            a(i,noutcom) = a(i,nft)
            a(i,nft) = atmp
        END DO
    END IF
    a4=a
    
    CALL cpu_time(start_time)
    
    IF (iscale == 1) CALL scaling

    WRITE(42,*)
    WRITE(42,*) 'Number of records: ',nrecord
    WRITE(42,*) 'Size of the training set: ',ntrain
    WRITE(42,*) 'Number of features: ',nft
    WRITE(41,*)
    WRITE(41,*) 'Number of records: ',nrecord
    WRITE(41,*) 'Size of the training set: ',ntrain
    WRITE(41,*) 'Number of features: ',nft
        
    x_0 = 2_dp  ! the very first x_0
    q=1         ! for linear regression
    maxdim=nft  ! for linear regression
    maxmin2 = 1 ! for linear regression
    mk=1        ! for linear regression 
    counter=0
    
    ! First compute just linear regression
    ! That is, we have only one linear function under min and only one min under max.
    ! The number of variables in the optimization problem is nft.
        
    CALL DBDC_algorithm( x_0(1:maxdim), x_solution(1:maxdim), f_solution, mit, mrounds, &
        & mrounds_clarke, termination, counter, agg_used,  &
        & stepsize_used, iprint )
        
    neval=neval+counter(3)+counter(7)    
    WRITE(42,*)
    WRITE(42,*) 'Fit function value for linear regression: ',f_solution
    WRITE(42,*) 'Number of linear function evaluations: ',neval
    IF (iscale == 1) THEN 
        CALL rescaling !(x_solution,x_unscaled)
    ELSE
        x_unscaled=x_solution
    END IF
    WRITE(42,*)

     
    ! Print the solution    
    WRITE(42,*) 'x: ',x_unscaled(1:maxdim) !
    x0reg=x_solution ! If scaling is used this is the scaled starting point.
    f_1=f_solution

    ! Indices for training data for linear regression

    IF (nrecord > ntrain .AND. ipredi == 1) THEN ! there are training and test sets
        e0_tr=0_dp
        ytest=0_dp
        DO i = 1, ntrain
            e0_tr=e0_tr+a(i,nft)
            pred(i)=x_unscaled(nft)
            DO j=1, nft-1
                pred(i) = pred(i)+x_unscaled(j)*a(i,j)
            END DO
            ytest=ytest+pred(i)
        END DO
        ytest=ytest/REAL(ntrain,dp)
        e0_tr=e0_tr/REAL(ntrain,dp)
        rmse1=0_dp
        mae1=0_dp
        e1_tr=0_dp
        DO i=1,ntrain
            rmse1=rmse1+(pred(i)-a(i,nft))**2
            mae1=mae1+ABS(pred(i)-a(i,nft))
            e1_tr=e1_tr+(a(i,nft)-e0_tr)**2 !e0?
        END DO
        mse1=rmse1/REAL(ntrain-nft,dp)
        ce1=1_dp-(rmse1/e1_tr)
        rmse1=SQRT(rmse1/REAL(ntrain,dp))
        mae1=mae1/REAL(ntrain,dp)


        IF(ntrain > 1) THEN 
            sx=0_dp
            sy=0_dp
            rcor=0_dp
            DO i=1,ntrain
                sx=sx+(pred(i)-ytest)**2
                sy=sy+(a(i,nft)-e0_tr)**2 
                rcor=rcor+(pred(i)-ytest)*(a(i,nft)-e0_tr)
            END DO
            sx=SQRT(sx/REAL(ntrain-1))
            sy=SQRT(sy/REAL(ntrain-1))
            
            r1=rcor/(REAL(ntrain-1,dp)*sx*sy)
        END IF

        WRITE(41,*)
        WRITE(41,*) 'Evaluation criteria for training set in linear regression (rmse, mae, mse, ce, r): '
        WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
        WRITE(41,*) mae1
        WRITE(41,*) mse1
        WRITE(41,*) ce1
        WRITE(41,*) r1


    END IF

    ! Indices for test data for linear regression
    
    IF (nrecord > ntrain .AND. ipredi == 1) THEN ! prediction
        ntest=nrecord-ntrain
        
        e0_test=0_dp
        ytest=0_dp
        DO i = 1, ntest
            e0_test=e0_test+a(ntrain+i,nft)
            pred(i)=x_unscaled(nft)
            DO j=1, nft-1
                pred(i) = pred(i)+x_unscaled(j)*a(ntrain+i,j)
            END DO    
            ytest=ytest+pred(i)
        END DO
        ytest=ytest/REAL(ntest,dp)
        e0_test=e0_test/REAL(ntest,dp)
        rmse1=0_dp
        mae1=0_dp
        e1_test=0_dp
        DO i=1,ntest
            rmse1=rmse1+(pred(i)-a(ntrain+i,nft))**2
            mae1=mae1+ABS(pred(i)-a(ntrain+i,nft))
            e1_test=e1_test+(a(ntrain+i,nft)-e0_test)**2
        END DO
        mse1=rmse1/REAL(ntest-nft,dp)
        ce1=1_dp-(rmse1/e1_test)
        rmse1=SQRT(rmse1/REAL(ntest,dp))
        mae1=mae1/REAL(ntest,dp)
        
        
        IF(ntest > 1) THEN
            sx=0_dp
            sy=0_dp
            rcor=0_dp
            DO i=1,ntest
                sx=sx+(pred(i)-ytest)**2
                sy=sy+(a(ntrain+i,nft)-e0_test)**2
                rcor=rcor+(pred(i)-ytest)*(a(ntrain+i,nft)-e0_test)
            END DO
            sx=SQRT(sx/REAL(ntest-1))
            sy=SQRT(sy/REAL(ntest-1))       
            
            r1=rcor/(REAL(ntest-1,dp)*sx*sy)
        END IF
        
        WRITE(41,*)
        WRITE(41,*) 'Evaluation criteria for linear regression (rmse, mae, mse, ce, r): ' 
        WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
        WRITE(41,*) mae1
        WRITE(41,*) mse1
        WRITE(41,*) ce1
        WRITE(41,*) r1
        
    
    END IF
        
    f_2=f_1
        
    IF (inc >= 1) THEN ! Compute starting points and structure incremental.
    
        ! Computation of min
    
        minloop: DO j=2,maxlin 
            PRINT*,j,'linear functions under min.'
            ! Set a new starting point
            x_0(1:maxdim)=x_solution(1:maxdim) ! maxdim first elements, note that maxdim is not yet updated.
            x_0(maxdim+1:j*nft)=x_solution(maxdim-nft+1:maxdim) ! the last nft elements
            mk=mk+1
            maxdim=maxdim+nft
            q=q+1
        
            CALL DBDC_algorithm( x_0(1:maxdim), x_solution(1:maxdim), f_solution, mit, mrounds, &
                & mrounds_clarke, termination, counter, agg_used,  &
                & stepsize_used, iprint )
        
            neval=neval + j*(counter(3)+counter(7))
            WRITE(42,*)
            WRITE(42,101) 1,j
            WRITE(42,*) 'Fit function value for PWLR is ',f_solution
            WRITE(42,*) 'Number of linear function evaluations: ',neval
            WRITE(42,*) 'The difference between function values is ',(f_2-f_solution)/(f_1+1_dp)
            WRITE(42,*) 'The 2. diff. between function values is   ',(f_2-f_solution)/(f_2+1_dp)
            WRITE(42,*)
101         FORMAT(' Model with maxmin = ',i3,' and maxlin = ', i3)
        
            
            
            IF (iscale == 1) THEN 
                CALL rescaling !(x_solution,x_unscaled)
            ELSE
                x_unscaled=x_solution
            END IF
            ! Print the solution    
            WRITE(42,*) 'x: ',x_unscaled(1:maxdim)
               
            ! Stop here if f is small enough   
            IF (f_solution < tol) EXIT minloop
            
            ! Commented, causes problems with next loop arrays.
            ! If you add this you should adjust maxlin etc corresbondingly.
            !IF ((f_2-f_solution)/(f_1+1_dp) < tol) THEN
            !PRINT*,' Early stopping from the first loop. ' 
            !EXIT minloop
            !END IF 
        
            f_2=f_solution
            
            
            ! Indices for training data for PWL-regression 

            IF (nrecord > ntrain .AND. ipredi == 1) THEN ! training and test sets are used
                WRITE(41,*)
                WRITE(41,101) 1,j
                mtmp=nft*j


                ytest=0_dp
                DO i = 1, ntrain
                    ! computation of correct piece
                    pind=1
                    ipiecebest=1
                    piece_best=-large
                    ipiecemin=pind
                    piece_min=large
                    DO j2=1,j
                        piece=x_unscaled(pind*nft) !
                        DO j3=1, nft-1
                            piece = piece+x_unscaled((pind-1)*nft+j3)*a(i,j3)
                        END DO

                        IF (piece<piece_min) THEN
                            ipiecemin=pind
                            piece_min=piece
                        END IF

                        pind=pind+1
                    END DO

                    IF (piece_min > piece_best) THEN
                        ipiecebest=ipiecemin
                        piece_best=piece_min
                    END IF

                    pred(i)=x_unscaled(ipiecebest*nft)
                    DO j1=1, nft-1
                        pred(i) = pred(i)+x_unscaled((ipiecebest-1)*nft+j1)*a(i,j1)
                    END DO
                    ytest=ytest+pred(i)
                END DO
                ytest=ytest/REAL(ntrain,dp)
                rmse1=0_dp
                mae1=0_dp
                DO i=1,ntrain
                    rmse1=rmse1+(pred(i)-a(i,nft))**2
                    mae1=mae1+ABS(pred(i)-a(i,nft))
                END DO
                mse1=rmse1/REAL(ntrain-mtmp,dp)
                ce1=1_dp-(rmse1/e1_tr)
                rmse1=SQRT(rmse1/REAL(ntrain,dp))
                mae1=mae1/REAL(ntrain,dp)

                IF(ntrain > 1) THEN
                    sx=0_dp
                    sy=0_dp
                    rcor=0_dp
                    DO i=1,ntrain
                        sx=sx+(pred(i)-ytest)**2
                        sy=sy+(a(i,nft)-e0_tr)**2 
                        rcor=rcor+(pred(i)-ytest)*(a(i,nft)-e0_tr)
                    END DO
                    sx=SQRT(sx/REAL(ntrain-1))
                    sy=SQRT(sy/REAL(ntrain-1))
                        
                    r1=rcor/(REAL(ntrain-1,dp)*sx*sy)
                END IF

                WRITE(41,*)
                WRITE(41,*) 'Evaluation criteria for training set in SVM piecewise linear regression (rmse, mae, mse, ce, r): '
                WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
                WRITE(41,*) mae1
                WRITE(41,*) mse1
                WRITE(41,*) ce1
                WRITE(41,*) r1
            END IF
            
            
            ! Indices for test data for PWL-regression 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF (nrecord > ntrain .AND. ipredi == 1) THEN ! prediction
                WRITE(41,*)
                WRITE(41,101) 1,j
                mtmp=nft*j
        
        
                ytest=0_dp
                DO i = 1, ntest
                    ! computation of correct piece
                    pind=1
                    ipiecebest=1
                    piece_best=-large
                    ipiecemin=pind
                    piece_min=large
                    DO j2=1,j
                        piece=x_unscaled(pind*nft) ! 
                        DO j3=1, nft-1
                            piece = piece+x_unscaled((pind-1)*nft+j3)*a(ntrain+i,j3)
                        END DO
                            
                        IF (piece<piece_min) THEN
                            ipiecemin=pind
                            piece_min=piece
                        END IF
                        
                        pind=pind+1
                    END DO 
                        
                    IF (piece_min > piece_best) THEN
                        ipiecebest=ipiecemin
                        piece_best=piece_min
                    END IF
                         
                  
                    pred(i)=x_unscaled(ipiecebest*nft)
                    DO j1=1, nft-1
                        pred(i) = pred(i)+x_unscaled((ipiecebest-1)*nft+j1)*a(ntrain+i,j1)
                    END DO    
                    ytest=ytest+pred(i)
                END DO
                ytest=ytest/REAL(ntest,dp)
                rmse1=0_dp
                mae1=0_dp
                DO i=1,ntest
                    rmse1=rmse1+(pred(i)-a(ntrain+i,nft))**2
                    mae1=mae1+ABS(pred(i)-a(ntrain+i,nft))
                END DO
                mse1=rmse1/REAL(ntest-mtmp,dp)
                ce1=1_dp-(rmse1/e1_test)
                rmse1=SQRT(rmse1/REAL(ntest,dp))
                mae1=mae1/REAL(ntest,dp)
                
                IF(ntest > 1) THEN
                    sx=0_dp
                    sy=0_dp
                    rcor=0_dp
                    DO i=1,ntest
                        sx=sx+(pred(i)-ytest)**2
                        sy=sy+(a(ntrain+i,nft)-e0_test)**2 
                        rcor=rcor+(pred(i)-ytest)*(a(ntrain+i,nft)-e0_test)
                    END DO
                    sx=SQRT(sx/REAL(ntest-1))
                    sy=SQRT(sy/REAL(ntest-1))       
                        
                    r1=rcor/(REAL(ntest-1,dp)*sx*sy)
                    
                END IF
                
                WRITE(41,*)
                WRITE(41,*) 'Evaluation criteria for SVM piecewise linear regression (rmse, mae, mse, ce, r): ' 
                WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
                WRITE(41,*) mae1
                WRITE(41,*) mse1
                WRITE(41,*) ce1
                WRITE(41,*) r1
            END IF
            

        END DO minloop
        
        
        IF (f_solution > tol) THEN
          
            maxloop: DO j=2, maxmin
                PRINT*,j,'min under max'
                f_2=f_solution 
            
                ! Set a new starting point
                x_0(1:maxdim)=x_solution(1:maxdim) ! maxdim first elements, note that maxdim is not yet updated.
                x_0(maxdim+1:j*maxlin*nft)=x_solution(maxdim-nft*maxlin+1:maxdim) ! the last nft elements
                
                
                q=maxlin*j
                maxdim=nft*q
                maxmin2=j
            
                
                CALL DBDC_algorithm( x_0(1:maxdim), x_solution(1:maxdim), f_solution, mit, mrounds, &
                    & mrounds_clarke, termination, counter, agg_used,  &
                    & stepsize_used, iprint )
                
                neval = neval+q*(counter(3)+counter(7))

                IF (inc==2) THEN ! Here we need to compute more if inc==2
                
                    ! Set a new starting point 
                    DO i=1,j*maxlin
                        x_0((i-1)*nft+1:i*nft)=x0reg(1:nft)
                    END DO

                    CALL DBDC_algorithm( x_0(1:maxdim), x_solution2(1:maxdim), f_solution2, mit, mrounds, &
                        & mrounds_clarke, termination, counter, agg_used,  &
                        & stepsize_used, iprint )

                    neval = neval+q*(counter(3)+counter(7))

                    IF (f_solution2 < f_solution) THEN
                        x_solution = x_solution2
                        f_solution = f_solution2
                    END IF

                END IF

                WRITE(42,*)
                WRITE(42,101) j,maxlin
                WRITE(42,*) 'Fit function value1 for PWLR is ',f_solution
                !WRITE(42,*) 'Fit function value2 for PWLR is ',f_solution2
                WRITE(42,*) 'Number of linear function evaluations: ',neval,q
                WRITE(42,*) 'The difference between function values is ',(f_2-f_solution)/(f_1+1_dp)
                !WRITE(42,*) 'The difference between function values is ',(f_2-f_solution2)/(f_1+1_dp)
                WRITE(42,*) 'The 2. diff. between function values is   ',(f_2-f_solution)/(f_2+1_dp)
                WRITE(42,*)

                !IF (inc==2) THEN ! This could be given above but then we cannot print. Fix to final version.

                !    IF (f_solution2 < f_solution) THEN
                !        x_solution = x_solution2
                !        f_solution = f_solution2
                !    END IF

                !END IF

            
                IF (iscale == 1) THEN 
                    CALL rescaling !(x_solution,x_unscaled)
                ELSE
                    x_unscaled=x_solution
                END IF
    
                !result can be printed here
                WRITE(42,*) 'x: ',x_unscaled(1:maxdim)
                
                
                IF (nft == 2 .AND. maxlin == 3) THEN ! Note works only for two dimensional data and 3 linear functions under maximum
                
                    WRITE(43,101) j,maxlin
                    h=minv
                    h2=(maxv-minv)/1000_dp
                    DO i=1,1000
                        tmp2 = - large
                        h=h+h2
                        i3=1
                        DO i2=1,j
                            tmp = MIN(x_unscaled(i3)*h+x_unscaled(i3+1),x_unscaled(i3+2)*h+x_unscaled(i3+3))
                            tmp = MIN(tmp,x_unscaled(i3+4)*h+x_unscaled(i3+5))
                            i3 = i3 + 6
                            IF (tmp > tmp2) tmp2 = tmp
                        END DO
                        WRITE(43,*) h, tmp2

                    END DO
                END IF

                mtmp=j  
                
                ! Stop here if no progress   
                
                !IF (f_solution < tol) EXIT maxloop 
                !IF ((f_2-f_solution)/(f_1+1_dp) < tol) EXIT maxloop 
        
                ! Indices for training

                IF (nrecord > ntrain .AND. ipredi == 1) THEN ! prediction
                    WRITE(41,*)
                    WRITE(41,101) mtmp,maxlin
                    mtmp=mtmp*nft*maxlin


                    ytest=0_dp
                    DO i = 1, ntrain
                        ! computation of correct piece
                        pind=1
                        ipiecebest=1
                        piece_best=-large
                        DO j1=1,j
                            ipiecemin=pind
                            piece_min=large
                            DO j2=1,maxlin
                                piece=x_unscaled(pind*nft) !
                                DO j3=1, nft-1
                                    piece = piece+x_unscaled((pind-1)*nft+j3)*a(i,j3)
                                END DO

                                IF (piece<piece_min) THEN
                                    ipiecemin=pind
                                    piece_min=piece
                                END IF

                                pind=pind+1
                            END DO

                            IF (piece_min > piece_best) THEN
                                ipiecebest=ipiecemin
                                piece_best=piece_min
                            END IF

                        END DO

                        pred(i)=x_unscaled(ipiecebest*nft)
                        DO j1=1, nft-1
                            pred(i) = pred(i)+x_unscaled((ipiecebest-1)*nft+j1)*a(i,j1)
                        END DO
                        ytest=ytest+pred(i)
                    END DO
                    ytest=ytest/REAL(ntrain,dp)
                    rmse1=0_dp
                    mae1=0_dp
                    DO i=1,ntrain
                        rmse1=rmse1+(pred(i)-a(i,nft))**2
                        mae1=mae1+ABS(pred(i)-a(i,nft))
                    END DO
                    mse1=rmse1/REAL(ntrain-mtmp,dp)
                    ce1=1_dp-(rmse1/e1_tr) !
                    rmse1=SQRT(rmse1/REAL(ntrain,dp))
                    mae1=mae1/REAL(ntrain,dp)

                    IF(ntrain > 1) THEN
                        sx=0_dp
                        sy=0_dp
                        rcor=0_dp
                        DO i=1,ntrain
                            sx=sx+(pred(i)-ytest)**2
                            sy=sy+(a(i,nft)-e0_tr)**2 
                            rcor=rcor+(pred(i)-ytest)*(a(i,nft)-e0_tr)
                        END DO
                        sx=SQRT(sx/REAL(ntrain-1))
                        sy=SQRT(sy/REAL(ntrain-1))
                        
                        r1=rcor/(REAL(ntrain-1,dp)*sx*sy)
                    END IF

                    WRITE(41,*)
                    WRITE(41,*) 'Evaluation criteria for training in SVM piecewise linear regression (rmse, mae, mse, ce, r): '
                    WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
                    WRITE(41,*) mae1
                    WRITE(41,*) mse1
                    WRITE(41,*) ce1
                    WRITE(41,*) r1
                END IF


                ! Indices for test set.

                mtmp=j
                IF (nrecord > ntrain .AND. ipredi == 1) THEN ! prediction
                    WRITE(41,*)
                    WRITE(41,101) mtmp,maxlin
                    mtmp=mtmp*nft*maxlin
        
        
                    ytest=0_dp
                    DO i = 1, ntest
                        ! computation of correct piece
                        pind=1
                        ipiecebest=1
                        piece_best=-large
                        DO j1=1,j
                            ipiecemin=pind
                            piece_min=large
                            DO j2=1,maxlin
                                piece=x_unscaled(pind*nft) ! 
                                DO j3=1, nft-1
                                    piece = piece+x_unscaled((pind-1)*nft+j3)*a(ntrain+i,j3)
                                END DO
                            
                                IF (piece<piece_min) THEN
                                    ipiecemin=pind
                                    piece_min=piece
                                END IF
                        
                                pind=pind+1
                            END DO 
                        
                            IF (piece_min > piece_best) THEN
                                ipiecebest=ipiecemin
                                piece_best=piece_min
                            END IF
                         
                        END DO
                
                        pred(i)=x_unscaled(ipiecebest*nft)
                        DO j1=1, nft-1
                            pred(i) = pred(i)+x_unscaled((ipiecebest-1)*nft+j1)*a(ntrain+i,j1)
                        END DO    
                        ytest=ytest+pred(i)
                    END DO
                    ytest=ytest/REAL(ntest,dp)
                    rmse1=0_dp
                    mae1=0_dp
                    DO i=1,ntest
                        rmse1=rmse1+(pred(i)-a(ntrain+i,nft))**2
                        mae1=mae1+ABS(pred(i)-a(ntrain+i,nft))
                    END DO
                    mse1=rmse1/REAL(ntest-mtmp,dp)
                    ce1=1_dp-(rmse1/e1_test)
                    rmse1=SQRT(rmse1/REAL(ntest,dp))
                    mae1=mae1/REAL(ntest,dp)
                
                    IF(ntest > 1) THEN
                        sx=0_dp
                        sy=0_dp
                        rcor=0_dp
                        DO i=1,ntest
                            sx=sx+(pred(i)-ytest)**2
                            sy=sy+(a(ntrain+i,nft)-e0_test)**2 
                            rcor=rcor+(pred(i)-ytest)*(a(ntrain+i,nft)-e0_test)
                        END DO
                        sx=SQRT(sx/REAL(ntest-1))
                        sy=SQRT(sy/REAL(ntest-1))       
                        
                        r1=rcor/(REAL(ntest-1,dp)*sx*sy)
                    
                    END IF
                
                    WRITE(41,*)
                    WRITE(41,*) 'Evaluation criteria for SVM piecewise linear regression (rmse, mae, mse, ce, r): ' 
                    WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
                    WRITE(41,*) mae1
                    WRITE(41,*) mse1
                    WRITE(41,*) ce1
                    WRITE(41,*) r1
                END IF
        
                mtmp=j
                IF (f_solution < tol) EXIT maxloop 
                IF ((f_2-f_solution)/(f_1+1_dp) < tol) EXIT maxloop 
        
                
            END DO maxloop
            
            IF (ipredi == 0) THEN ! Print only the final result
                IF (nrecord > ntrain) THEN ! prediction
                    WRITE(41,*)
                    WRITE(41,*) 'Final model with maxmin = ',mtmp,' and maxlin = ', maxlin
                    mtmp=mtmp*nft*maxlin
        
        
                    ytest=0_dp
                    DO i = 1, ntest
                        ! computation of correct piece
                        pind=1
                        ipiecebest=1
                        piece_best=-large
                        DO j=1,maxmin
                            ipiecemin=pind
                            piece_min=large
                            DO j2=1,maxlin
                                piece=x_unscaled(pind*nft) ! 
                                DO j3=1, nft-1
                                    piece = piece+x_unscaled((pind-1)*nft+j3)*a(ntrain+i,j3)
                                END DO
                            
                                IF (piece<piece_min) THEN
                                    ipiecemin=pind
                                    piece_min=piece
                                END IF
                        
                                pind=pind+1
                            END DO 
                        
                            IF (piece_min > piece_best) THEN
                                ipiecebest=ipiecemin
                                piece_best=piece_min
                            END IF
                         
                        END DO
                
                        pred(i)=x_unscaled(ipiecebest*nft)
                        DO j=1, nft-1
                            pred(i) = pred(i)+x_unscaled((ipiecebest-1)*nft+j)*a(ntrain+i,j)
                        END DO    
                        ytest=ytest+pred(i)
                    END DO
                    ytest=ytest/REAL(ntest,dp)
                    rmse1=0_dp
                    mae1=0_dp
                    DO i=1,ntest
                        rmse1=rmse1+(pred(i)-a(ntrain+i,nft))**2
                        mae1=mae1+ABS(pred(i)-a(ntrain+i,nft))
                    END DO
                    mse1=rmse1/REAL(ntest-mtmp,dp)
                    ce1=1_dp-(rmse1/e1_test)
                    rmse1=SQRT(rmse1/REAL(ntest,dp))
                    mae1=mae1/REAL(ntest,dp)
                
                    IF(ntest > 1) THEN
                        sx=0_dp
                        sy=0_dp
                        rcor=0_dp
                        DO i=1,ntest
                            sx=sx+(pred(i)-ytest)**2
                            sy=sy+(a(ntrain+i,nft)-e0_test)**2 
                            rcor=rcor+(pred(i)-ytest)*(a(ntrain+i,nft)-e0_test)
                        END DO
                        sx=SQRT(sx/REAL(ntest-1))
                        sy=SQRT(sy/REAL(ntest-1))       
                        
                        r1=rcor/(REAL(ntest-1,dp)*sx*sy)
                   
                    END IF
                
                    WRITE(41,*)
                    WRITE(41,*) 'Evaluation criteria for SVM piecewise linear regression (rmse, mae, mse, ce, r): ' 
                    WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
                    WRITE(41,*) mae1
                    WRITE(41,*) mse1
                    WRITE(41,*) ce1
                    WRITE(41,*) r1
        
                END IF
            END IF
            
        
        
        ELSE ! f_solution <= tol with only on min function under the maximum
            WRITE(42,*)
            WRITE(42,*) 'Fit function value for piecewise linear regression with ', j
            WRITE(42,*) 'linear functions under the minimum and ',1,' min under the maximum'
            WRITE(42,*) 'is ',f_solution
            IF (iscale == 1) THEN 
                CALL rescaling !(x_solution,x_unscaled)
            ELSE
                x_unscaled=x_solution
            END IF
            
            ! Print the solution    
            WRITE(42,*) 'x: ',x_unscaled
               
        
            WRITE(41,*)
            WRITE(41,*) 'Final model with maxmin = ',1,' and maxlin = ', j
            WRITE(41,*) 'Write prediction to "else" part.'
            WRITE(41,*) 'This never happened in our numerical experiments.' 
            WRITE(41,*) 'If you get this result, please, send message to napsu@karmitsa.fi'
        END IF
        
        
    ELSE ! Compute the known model with maxmin and maxlin functions
        ! Set a new starting point
        DO i=1,maxmin*maxlin 
            x_0((i-1)*maxdim+1:i*maxdim)=x_solution(1:maxdim)
        END DO
        
        CALL init_par() ! sets dimension for optimization problem and q.
        mk=maxlin ! each min function has maxlin linear functions
        maxmin2=maxmin
             
        CALL DBDC_algorithm( x_0, x_solution, f_solution, mit, mrounds, &
            & mrounds_clarke, termination, counter, agg_used,  &
            & stepsize_used, iprint )
        
        WRITE(42,*)
        WRITE(42,*) 'Fit function value for PWLR: ',f_solution
        IF (iscale == 1) THEN 
            CALL rescaling !(x_solution,x_unscaled)
        ELSE
            x_unscaled=x_solution
        END IF
    
        ! print the results
        WRITE(42,*) 'x: ',x_unscaled
        
        IF (nrecord > ntrain) THEN ! prediction
            WRITE(41,*)
            WRITE(41,*) 'Final model with maxmin = ',maxmin,' and maxlin = ', maxlin
            
            ytest=0_dp
            DO i = 1, ntest
                ! computation of correct piece
                pind=1
                ipiecebest=1
                piece_best=-large
                DO j=1,maxmin
                    ipiecemin=pind
                    piece_min=large
                    DO j2=1,maxlin
                        piece=x_unscaled(pind*nft) ! 
                        DO j3=1, nft-1
                            piece = piece+x_unscaled((pind-1)*nft+j3)*a(ntrain+i,j3)
                        END DO
                            
                        IF (piece<piece_min) THEN
                            ipiecemin=pind
                            piece_min=piece
                        END IF
                        
                        pind=pind+1
                    END DO 
                        
                    IF (piece_min > piece_best) THEN
                        ipiecebest=ipiecemin
                        piece_best=piece_min
                    END IF
                         
                END DO
                
                pred(i)=x_unscaled(ipiecebest*nft)
                DO j=1, nft-1
                    pred(i) = pred(i)+x_unscaled((ipiecebest-1)*nft+j)*a(ntrain+i,j)
                END DO   
                ytest=ytest+pred(i) 
            END DO
            ytest=ytest/REAL(ntest,dp)
            rmse1=0_dp
            mae1=0_dp
            DO i=1,ntest
                rmse1=rmse1+(pred(i)-a(ntrain+i,nft))**2
                mae1=mae1+ABS(pred(i)-a(ntrain+i,nft))
            END DO
            mse1=rmse1/REAL(ntest-nft*maxmin*maxlin,dp)
            ce1=1_dp-(rmse1/e1_test)
            rmse1=SQRT(rmse1/REAL(ntest,dp))
            mae1=mae1/REAL(ntest,dp)
            
            IF(ntest > 1) THEN
                sx=0_dp
                sy=0_dp
                rcor=0_dp
                DO i=1,ntest
                    sx=sx+(pred(i)-ytest)**2
                    sy=sy+(a(ntrain+i,nft)-e0_test)**2 
                    rcor=rcor+(pred(i)-ytest)*(a(ntrain+i,nft)-e0_test)
                END DO
                sx=SQRT(sx/REAL(ntest-1))
                sy=SQRT(sy/REAL(ntest-1))       
                       
       
                r1=rcor/(REAL(ntest-1,dp)*sx*sy)
            END IF
       
            WRITE(41,*)
            WRITE(41,*) 'Evaluation criteria for SVM piecewise linear regression (rmse, mae, mse, ce, r): ' 
            WRITE(41,*) rmse1!,mae1,mse1,ce1,r1
            WRITE(41,*) mae1
            WRITE(41,*) mse1
            WRITE(41,*) ce1
            WRITE(41,*) r1
        
        END IF
          
        
    END IF
        
    
    CALL cpu_time(finish_time)
    WRITE(41,*)
    WRITE(41,*) 'CPU time used: ',finish_time-start_time 
    WRITE(42,*)
    WRITE(42,*) 'CPU time used: ',finish_time-start_time 
           
    
    ! Changing the order back
    ! Maybe no need to do this. Instead do prediction.
    IF (noutcom < nft) THEN

        DO i=1,nrecord
            atmp = a(i,noutcom)
            a(i,noutcom)=a(i,nft)
            a(i,nft)=atmp
        END DO
    END IF

    

    CLOSE(42)
    CLOSE(41)
    CLOSE(43)

    CLOSE(78)
    deallocate(x_0,x0reg,x_solution,x_unscaled)
    deallocate(a,a4)
    deallocate(pred,dlin,clin)
    deallocate(mk)
    

    STOP


END PROGRAM spr
