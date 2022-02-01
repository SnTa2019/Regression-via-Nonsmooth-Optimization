!*************************************************************************
!*                                                                       *
!*     Initialization of parameters for A-PWLSVR                         *
!*     (last modified 01.02.2022)                                        *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     initspr          ! Initialization of parameters for A-PWLSVR (main program svr).
!*     initdbdc         ! Initialization of DBDC -solver.
!*



MODULE initspr  ! Initialization of parameters.

    USE constants, ONLY : dp   ! Precision for reals.
    IMPLICIT NONE

    ! Names of input and output files:
    CHARACTER(LEN=80), SAVE :: &
        infile, &                        ! Name of the dataset file. Given as commandline argument.
        outfile1 = 'reg_results.txt', &  ! Regression results with hyperplanes.
        outfile2 = 'predictions.txt', &  ! Regression results with predictions.
        outfile3 = 'pwl_points.txt'      ! Points of original data and pwl-model obtained
                                         ! for drawing purposis. Not working with nft > 2. 
        
    ! Input real parameters. 
    REAL(KIND=dp), SAVE :: &
        minv = -2_dp, &                  ! Minimum value for pwl_points (need to be given by the user), 
                                         ! used with outfile3, othewise ignored.
        maxv = 2_dp, &                   ! Maximum value for pwl_points (need to be given by the user), 
                                         ! used with outfile3, othewise ignored.
        tol = 0.05_dp                    ! Stopping tolerance.
        
    ! Input integer parameters.
    INTEGER, SAVE :: &
        nft, &                           ! Number of features in data, from user.
        nrecord, &                       ! Number of records in data, from user.
        ntrain, &                        ! Size of the training set, ntrain <= nrecord.
        noutcom, &                       ! Output column.
        maxmin = 2, &                    ! Maximum number of min functions in the model, from user.
        maxlin = 3, &                    ! Number of linear functions in the model, optional, from user.
        iscale = 1, &                    ! Scaling is
                                         !   0 - not used,
                                         !   1 - used (recommended).
        inc = 2, &                       ! Incremental algorithm is
                                         !   0 - not used (the starting point is obtained by regression),
                                         !   1 - used.
                                         !   2 - used with double starting points.
        ipredi = 1                       ! Prediction:
                                         !   0 - only at the termination with complete model.
                                         !   1 - every step.


    ! Allocatable real arrays.
    REAL(KIND=dp), SAVE, DIMENSION(:,:), allocatable :: &
        a, &                             ! Data set, from input file.
        a4                               ! Training part of the data (ntrain first objects of a).
        
      
    REAL(KIND=dp), SAVE, DIMENSION(:), allocatable :: & !
        x0, &                            ! the starting point. 
        x0reg, &                         ! SVR solution.
        x_solution, &                    ! the solution.
        x_solution2, &                   ! the solution.
        x_unscaled, &                    ! the unscaled solution.
        pred, &                          ! the values of the predicted outputs
        clin, &
        dlin
                                        
                                        
    ! Other real parameters.
    REAL(KIND=dp), SAVE :: &
        f_solution, &                    ! the objective function value at the solution 'x_solution'       
        f_solution2, &                   ! the objective function value at the solution 'x_solution2'
        eps = 0.050_dp, &                ! the epsilon used in SVM, change in init_par.
        CC = 0.50_dp, &                  ! the coefficient in front of maximum terms  (i.e. maximum from minimums)
        Cnull = 0.0_dp                   ! the coefficient in front of the quadratic terms, Should be zero!
        
    ! Other integer parameters.
    INTEGER, SAVE :: & !
        maxmin2, &                       ! maxmin2 = 1,...,maxmin
        maxdim, &                        ! Maximum dimension in optimization, nft*maxmin*maxlin
        q                                ! Total number of linear pieces, maxmin*maxlin
        
    ! Allocatable integer tables
    INTEGER, SAVE, DIMENSION(:), allocatable :: &
        mk                               ! mk(maxmin) = number of lines in each minimum function
        
CONTAINS

    SUBROUTINE init_par()  ! Initialization of parameters.

        IMPLICIT NONE
        
        q = maxmin*maxlin                ! this if maxlin is constant for all k=1,...,maxmin
        maxdim = nft*q
        IF (iscale == 1) eps=0.05_dp
        
    END SUBROUTINE init_par

    !===================================================================
    !  Scaling
    !===================================================================

    SUBROUTINE scaling

        IMPLICIT NONE

        ! Local variables
        REAL(KIND=dp) ::      &
            dm,var
        INTEGER ::  &
            i,j

        ! Intrinsic Functions
        INTRINSIC ABS,SQRT

        DO i=1,nft
            dm=0_dp
            DO j=1,ntrain   
                dm=dm+a4(j,i)
            END DO
            dm=dm/REAL(ntrain,dp)
            var=0_dp
            DO j=1,ntrain
                var=var+(a4(j,i)-dm)**2
            END DO
            var=var/REAL(ntrain-1,dp)
            var=SQRT(var)

            IF (var >= 1.0E-05_dp) THEN
                clin(i)=1_dp/var
                dlin(i)=-dm/var
                DO j=1,ntrain
                    a4(j,i)=clin(i)*a4(j,i)+dlin(i)
                END DO
            ELSE
                IF (ABS(dm) <= 1_dp) THEN
                    clin(i)=1_dp
                    dlin(i)=0_dp
                ELSE
                    clin(i)=1_dp/dm
                    dlin(i)=0_dp
                END IF
                DO j=1,ntrain
                    a4(j,i)=clin(i)*a4(j,i)+dlin(i)
                END DO
            END IF
        END DO

        RETURN

    END SUBROUTINE scaling

    !===================================================================
    !  Rescaling
    !===================================================================

    SUBROUTINE rescaling 

        IMPLICIT NONE

        ! Local variables
        REAL(KIND=dp) ::      &
            d0,rabs
        INTEGER ::  &
            i,j,j1,k

        ! Intrinsic Functions
        INTRINSIC ABS

        rabs=ABS(clin(nft))
        IF (rabs > 1.0E-08_dp) THEN
            DO k=1,maxmin2 
                DO j=1,mk(k) ! This is correct as we always do upto maxlin linear functions
                    i=j+(k-1)*mk(k)
                !DO j=1,maxlin ! This is correct as we always do upto maxlin linear functions
                !    i=j+(k-1)*maxlin
                    d0=0_dp
                    DO j1=1,nft-1
                        x_unscaled(j1+(i-1)*nft)=x_solution(j1+(i-1)*nft)*clin(j1)/rabs
                        d0=d0+x_solution(j1+(i-1)*nft)*dlin(j1)
                    END DO
                    x_unscaled(i*nft)=(x_solution(i*nft)+d0-dlin(nft))/rabs
                END DO
            END DO
        END IF

        RETURN

    END SUBROUTINE rescaling

END MODULE initspr

MODULE initdbdc  ! Initialization of parameters for DBDC. There is, usually, no need to change these parameters.

    USE constants, ONLY : dp        ! Precision for reals.

    IMPLICIT NONE

    ! Parameters
    REAL(KIND=dp), PARAMETER :: &  
        user_m = 0.01_dp, &         ! The descent parameter. If user_m <= 0.0_dp .OR. user_m >= 1.0_dp
                                    ! then DEFAULT value 0.2_dp is used
        user_c = -0.99_dp, &        ! The decrease parameter c in DBDC. If user_c <= 0.0_dp or user_c > 1.0_dp   
                                    ! then DEFAULT value 0.1_dp is used 
        user_r_dec = -0.99_dp, &    ! The decrease parameter r in DBDC       
        user_r_inc = (10.0_dp)**(7), & ! The increase parameter R: If user_r_inc <= 1.0_dp 
                                    ! then DEFAULT value (10.0_dp)**7 is used
        user_eps_1 = 5*(10.0_dp)**(-5), & ! The enlargement parameter: If user_eps_1 <= 0.0_dp .OR. user_eps_1 > 1.0_dp 
                                    ! then DEFAULT value 5*(10.0_dp)**(-5) is used
        user_crit_tol = (10.0_dp)**(-3), & ! If user_crit_tol <= 0.0_dp then DEFAULT value (10.0_dp)**(-5) is used when n <=200
                                            ! (10.0_dp)**(-4) is used when n > 200
        user_m_clarke = 0.01_dp, &  ! The descent parameter: If user_m_clarke <= 0.0_dp .OR. user_m_clarke >= 1.0_dp 
                                    ! then DEFAULT value 0.01_dp is used
        user_eps = 10.0_dp**(-4)    ! The proximity measure: If user_eps <= 0.0_dp
                                    ! then DEFAULT value (10.0_dp)**(-6) is used when n <= 50
                                    ! (10.0_dp)**(-5) is used when n > 50
                                                                    
                                     
    INTEGER, PARAMETER :: &
        problem1 = 1, &             ! Objective function selected for the DC component f1
        problem2 = problem1, &      ! Objective function selected for the DC component f2
        user_size_b2 = 3            ! The biggest possible size of the bundle B_2
                                    ! If user_size_b2 <= 0 then DEFAULT value 3 is used  
       
    INTEGER, PARAMETER :: &
        user_n = 2, & 
        user_size_b1 = MIN(user_n + 3,50), & ! The biggest possible size of the bundle B_1 
                                    ! If user_size_b1 <= 0 then DEFAULT value MIN(user_n+5,1000) is used 
        user_size = user_n * 2      ! The biggest possible bundle size for the bundle used in 
                                    ! 'Clarke stationary' algorithm                                                            
  
    INTEGER, SAVE, DIMENSION(8) :: counter  ! Contains the values of different counteres:
                                    !   counter(1) = iter_counter:         the number of 'main iterations' executed
                                    !   counter(2) = subprob_counter:      the number of subproblems solved in 'main iteration'
                                    !   counter(3) = f_counter:            the number of function values evaluated for DC component in 'main iteration'
                                    !   counter(4) = subgrad1_counter:     the number of subgradients calculated for f_1 in 'main iteration'
                                    !   counter(5) = subgrad2_counter:     the number of subgradients calculated for f_2 in 'main iteration'
                                    !   counter(6) = stop_cond_counter:    the number of times 'Clarke stationary algorithm' is executed 
                                    !   counter(7) = clarke_f_counter:     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                    !   counter(8) = clarke_sub_counter:   the number of subgradients caluculated for f in 'Clarke stationary algorithms'
 
    
    INTEGER :: &
        iprint = 0, &               ! Variable that specifies print option (specified by USER): 
                                    !   iprint = 0 : print is suppressed
                                    !   iprint = 1 : basic print of final result 
                                    !   iprint = -1: basic print of final result (without the solution vector)
                                    !   iprint = 2 : extended print of final result 
                                    !   iprint = -2: extended print of final result (without the solution vector)
                                    !   iprint = 3 : basic print of intermediate results and extended print of final results
                                    !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                    !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                    !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                    !
                                    ! If 'iprint' <= -5 .OR. 'iprint' >= 5 then DEFAULT value 'iprint'=1 is used    
        mit = 500, &                ! The maximum number of 'main iterations' (specified by USER).
                                    ! If 'mit' <=0 then DEFAULT value 'mit'=1000 is used
        mrounds = 500, &            ! The maximum number of rounds during one 'main iteration' (specified by USER).
                                    ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used
        mrounds_clarke = 500, &     ! The maximum number of rounds during one 'Clarke stationarity' algorithm (specified by USER).
                                    ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
        termination                 ! The reason for termination in DBDC:
                                    !   1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                    !   2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < user_eps)
                                    !   3 - the maximum number 'mrounds' of rounds is executed in one main iteration
                                    !   4 - the maximum number of 'main iterations' is executed 
                                    !   5 - the maximum number 'mrounds_clarke' of rounds is executed in one 'Clarke stationary' algorithm

    LOGICAL :: &
        !agg_used = .FALSE., &      ! If .TRUE. then aggregation is used in DBDC (specified by USER).
        agg_used = .TRUE., &        ! If .TRUE. then aggregation is used in DBDC (specified by USER).
        stepsize_used = .TRUE.      ! If .TRUE. then simple stepsize determination is done in DBDC after each 'main iteration' (specified by USER).
                     
END MODULE initdbdc
