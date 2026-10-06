c=============================================================
c Incremental DGM algorithm for weighted clusterwise linear regression
c Hyperplanes found incrementally adding a hyperplane at each iteration
c Discrete gradient method is applied to solve optimization problems
c=============================================================
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Define variables, parameters, functions, ... used in the code

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! <><><><><><>< parameters ><><><><><>
! nft             :Number of features (including output/label)
! mft = nft-1 :Number of features without output
! nrecord      :Number of observations: 
! nc              :Number of clusters, hyperlines
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! <><><><><><><  variables ><><><><><>
! x

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! <><><><><><>< functions/subroutines ><><><><><>

! scaling          : scaling data including output
! rescaling       : reascaling data  
! clusterfunc   :Weighted CLR function 
! auxfunc        :Weighted Auxiliary CLR function 
! func 
! funcb

! combdecomp

! distribution
! checkfinal
! step1
! cleaning
! fgrad 
! optim 

! optimum
! wolfe 
! equations
! dgrad  
! dgrad2  
! armijo 
 
! gradient 
! bfgs 
! optimumb 
! dir 
! matrix1
! armijob

c===================================================================
c     main programm
c===================================================================
      PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100
     1 ,maxinit=10000, maxrec1=5000)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      DOUBLE PRECISION x(maxvar),a(maxrec,maxnft),a1(maxnft),x4(maxvar)
     1 ,dminim(maxrec),xc(maxinit,maxnft),z(maxvar),fv(maxrec)
     2 ,xbest(maxvar),a4(maxrec,maxnft),zc(maxinit,maxnft),z1(maxvar)
     7 ,zbest(maxvar),x1(maxvar),x2(maxvar),z2(maxvar),z4(maxvar)
     8 ,weight(maxnft),fc(maxclust),dist(maxrec1,maxrec1)
     3 ,ftr(maxclust,100),rmsetr(maxclust,100),rmsetest1(maxclust,100)
     4 ,rmsetest2(maxclust,100),rmsetest3(maxclust,100)
     5 ,rmsetest4(maxclust,100),rmsetest5(maxclust,100)
     5 ,rmsetest6(maxclust,100),rmaetr(maxclust,100)
     6 ,rmaetest1(maxclust,100),rmaetest2(maxclust,100)
     7 ,rmaetest3(maxclust,100),rmaetest4(maxclust,100)
     8 ,rmaetest5(maxclust,100),rmaetest6(maxclust,100)    
     9 ,cetr(maxclust,100), cetest1(maxclust,100),cetest2(maxclust,100)
     9 ,cetest3(maxclust,100),cetest4(maxclust,100)
     9 ,cetest5(maxclust,100),cetest6(maxclust,100),cortr(maxclust,100)
     9 ,cortest1(maxclust,100),cortest2(maxclust,100)
     9 ,cortest3(maxclust,100),cortest4(maxclust,100)
     9 ,cortest5(maxclust,100),cortest6(maxclust,100) 
     9 ,xx(maxnft),zz(maxnft)    
      INTEGER nob(maxclust,maxrec),nel(maxclust),nob1(maxrec)
     1 ,nlist(maxrec),listtest(maxrec,50)
      COMMON /c22/a,/anclust/nclust,/cnft/nft,/crecord/nrecord
     1 ,/cnk/nel,nob,/ctoler/toler,/cns/ns,/cgamma/gamma1,/cmft/mft
     2 ,/cdminim/dminim,/cnc/nc,/cnel1/nel1,nob1,/cnlist/nl,nlist
     3 ,/crecord1/nrecord1,/c24/a4,/cnrecord2/nrecord2,/ctlin/tlin
     4 ,/cnc1/nc1,/ccoef/coef,/cf1f2/f01,f02,/cweight/weight
     4 ,/cltest/knn,listtest,/cstab/tstab,/cdelta/delta
      OPEN(39,file='Results-inschgs.txt')
      OPEN(40,file='Statistics-inschgs.txt')
      OPEN(42,file='Function-value-inschgs.txt')
      OPEN(41,file='Hyperplane-inschgs.txt')
      OPEN(43,file='Predictions-inschgs.txt')
      OPEN(78,file='inschgs.txt',status='old',form='formatted')
c======================================================================
c               #attributes(input+output)       training
c  inschgs              6                         1071 
c    tstab=2.5;  balance=300
c======================================================================
      PRINT *,' '
      PRINT *,'Number of features (including output column):'
      READ *,nft
      PRINT *,' '
      PRINT *,'Output column:'
      READ *,noutcom
      PRINT *,' '
      PRINT *,'Number of clusters:'
      READ *,nclust
      PRINT *,' '
      PRINT *,'Value of stability parameter:'
      READ *,tstab
      PRINT *,' '
      PRINT *,'Value of balancing parameter:'
      READ *,balance
      PRINT *,' '
c====================================================================
      mft = nft - 1
      DO i=1,maxrec
       READ(78,*,END=901,ERR=900) (a(i,k),k=1,nft)
       nrecord=i
      END DO
  900 STOP 'Error in input file'       
  901 WRITE(40,*) 'Input complete. Number of records: ',nrecord
      WRITE(40,*)
c====================================================================
c      PRINT *,'Choose:'
c      PRINT *,'1: no test set'
      PRINT *,'2: one training set and one test set'
c      PRINT *,'3: cross-validation'
c      READ *,itest
      itest=2
c      IF(itest.eq.1) THEN
c        ntrain=nrecord
c        nfoldmax=1
c      END IF
c      IF(itest.eq.2) THEN
        PRINT *,'Enter the number of points in training set'
        READ *,ntrain
        nfoldmax=1
c      END IF
c      IF(itest.eq.3) THEN
c        PRINT *,'Enter the number of cross validations:'
c        READ *,nfoldmax
c      END IF
c======================================================================
      tlimit=7.2d+04
      delta=1.35d+00
      CALL cpu_time(time1)
    
      IF(noutcom.lt.nft) THEN  
       DO i=1,nrecord
        j1=0
        DO j=1,nft
         IF(j.ne.noutcom) THEN
          j1=j1+1
          a1(j1)=a(i,j)
         END IF
        END DO
        a1(nft)=a(i,noutcom)
        DO j=1,nft
         a(i,j)=a1(j)
        END DO
       END DO
      END IF
c====================================================================      
      IF(nrecord.le.200) THEN
       gamma1=0.0d+00
       gamma2=2.0d+00
      END IF

      IF((nrecord.gt.200).and.(nrecord.le.1000)) THEN
       gamma1=9.2d-01
       gamma2=1.2d+00
      END IF
      
      IF((nrecord.gt.1000).and.(nrecord.le.5000)) THEN
       gamma1=9.75d-01
       gamma2=1.02d+00
      END IF

      IF((nrecord.gt.5000).and.(nrecord.le.15000)) THEN
       gamma1=9.85d-01
       gamma2=1.02d+00
      END IF
      
      IF((nrecord.gt.15000).and.(nrecord.le.50000)) THEN
       gamma1=9.95d-01
       gamma2=1.02d+00
      END IF
      
      IF(nrecord.gt.50000) THEN
       gamma1=9.999d-01
       gamma2=1.002d+00
      END IF
c=================================================================
c Division to training and test sets (not ranDOmly)
c==================================================================
       write(42,*) 
       write(42,55) tstab
 55    format('The value of stability parameter:',f10.4)
       write(42,*) 
       write(42,56) balance
 56    format('The value of balancing parameter:',f10.4)
       write(42,*)
       write(42,*)
       write(42,*)
       WRITE(42,701)
 701   FORMAT('#Clust','             fvalue','        #Regprob_eval',
     1 '           #Lin_func_val','     #Lin_eval_per_point',
     2 '          CPU')
       WRITE(42,*)
c===================================================================
      nfold=0
  211 nfold=nfold+1
      IF(nfold.gt.nfoldmax) GO TO 210
      IF(itest.eq.1) THEN
       nl=0
       DO i=1,nrecord
        DO k=1,nft
         a4(i,k)=a(i,k)
        END DO
       END DO 
      END IF
      IF(itest.eq.2) THEN
       DO i=1,ntrain
        DO k=1,nft
         a4(i,k)=a(i,k)
        END DO
       END DO 
       nl=nrecord-ntrain
       DO i=1,nl
        nlist(i)=i+ntrain
       END DO
      END IF
      IF(itest.eq.3) THEN
       nrecord2=nrecord/nfoldmax
       nrecord3=(nfold-1)*nrecord2+1
       nrecord4=nfold*nrecord2
       nl=0
       DO j=nrecord3,nrecord4
        nl=nl+1
        nlist(j-nrecord3+1)=j 
       END DO
       nl1=0
       DO i=1,nrecord
        DO j=1,nl
         IF(i.eq.nlist(j)) GO TO 203
        END DO
        nl1=nl1+1
        DO k=1,nft
         a4(nl1,k)=a(i,k)
        END DO
 203   END DO
      END IF
      nrecord2=nl
      nrecord1=nrecord-nrecord2
c====================================================================      
      knn=nrecord1/100
      knn=max(knn,3)
      knn=min(knn,11)
      if((nl.gt.0).and.(nrecord.le.maxrec1)) then
       do i=1,nrecord
        dist(i,i)=0.0d+00
        do j=i+1,nrecord
         d3=0.0d+00
         do k=1,nft-1
          d3=d3+dabs(a(i,k)-a(j,k))
         end do
         dist(i,j)=d3
         dist(j,i)=dist(i,j)        
        end do
       end do      
c====================================================================
       do i=1,nrecord
        klist=0
 81     d6=1.0d+30
        do j=1,nrecord1
         if(j.eq.i) go to 51         
         do i2=1,klist
          IF(j.eq.listtest(i,i2)) go to 51
         end do
         d3=dist(i,j)
         if(d3.lt.d6) then
          d6=d3
          k1=j
         end if 
  51    end do
        klist=klist+1
        listtest(i,klist)=k1
        if(klist.lt.knn) go to 81
       end do
      end if
c--------------------------------------      
      if((nl.gt.0).and.(nrecord.gt.maxrec1)) then
       do i=1,nrecord
        klist=0
 181    d6=1.0d+30
        do j=1,nrecord1
         if(j.eq.i) go to 151
         do i2=1,klist
          IF(j.eq.listtest(i,i2)) go to 151
         end do
         d3=0.0d+00
         do k=1,nft-1
          d3=d3+dabs(a(i,k)-a(j,k))
         end do
         if(d3.lt.d6) then
          d6=d3
          k1=j
         end if 
 151    end do
        klist=klist+1
        listtest(i,klist)=k1
        if(klist.lt.knn) go to 181
       end do
      end if      
c====================================================================
      CALL scaling
c====================================================================
      ncount=0
      tlin=0.0d+00
      CALL cpu_time(time2)
      DO nc=1,nclust
       PRINT 3,nc
 3     FORMAT('         The number of hyperplanes:',i5)
 ! initialise variables x for CLR, z for clustering
       IF(nc.eq.1) THEN
        do j=1,mft
         weight(j)=1.0d+00/dble(mft)
        end do
        DO j=1,nft
         x(j)=1.0d+00
        END DO
        nel1=nrecord1
        DO k=1,nrecord1
         nob1(k)=k
        END DO
        do j=1,mft
         z(j)=0.0d+00
         do k=1,nrecord1
          z(j)=z(j)+a4(k,j)
         end do
         z(j)=z(j)/dble(nrecord1)
        end do
        CALL bfgs(x)
        do i=1,nft
         xx(i)=x(i)
        end do
        do i=1,mft
         zz(i)=z(i)
        end do
        coef=1.0d+00
        ncount=ncount+1
        GO TO 2
       END IF
       CALL step1(x,xc,zc,ngood0)                                       ! start with one linear function
       print *,ngood0
c=====================================================================
       fb=1.0d+32
       DO i=1,ngood0
        ns=1
        DO j=1,nft
         x1(j)=xc(i,j)
        END DO
        DO j=1,mft
         z1(j)=zc(i,j)
        END DO
        CALL optim(x1,z1,f)
        DO j=1,nft
         xc(i,j)=x1(j)
        END DO
        DO j=1,mft
         zc(i,j)=z1(j)
        END DO
        fv(i)=f
        fb=dmin1(fb,f)
       END DO

       fbest1=gamma2*fb
       jpoints=0
       DO i=1,ngood0
        IF(fv(i).le.fbest1) THEN
         jpoints=jpoints+1
         DO j=1,nft
          xc(jpoints,j)=xc(i,j)
         END DO
         DO j=1,mft
          zc(jpoints,j)=zc(i,j)
         END DO
        END IF
       END DO
       CALL cleaning(jpoints,xc,zc)
       ngood0=jpoints
c auxiliary problem finishes here
c222    continue
       fb=1.0d+32
       DO i=1,ngood0
        ns=2
        DO k=1,nc-1
         do j=1,nft
          x2(j+(k-1)*nft)=x(j+(k-1)*nft)
         end do
         do j=1,mft
          z2(j+(k-1)*nft)=z(j+(k-1)*nft)
         end do
        END DO
       
        DO j=1,nft
         x2(j+(nc-1)*nft)=xc(i,j)
        END DO
        DO j=1,mft
         z2(j+(nc-1)*mft)=zc(i,j)
        END DO
        CALL optim(x2,z2,f)
        IF(f.lt.fb) THEN
         fb=f
         DO k=1,nc
          DO j=1,nft
           xbest(j+(k-1)*nft)=x2(j+(k-1)*nft)
          END DO
          DO j=1,mft
           zbest(j+(k-1)*mft)=z2(j+(k-1)*mft)
          END DO
         END DO
        END IF
       END DO

       DO k=1,nc
        DO j=1,nft
         x(j+(k-1)*nft)=xbest(j+(k-1)*nft)
        END DO
        DO j=1,mft
         z(j+(k-1)*mft)=zbest(j+(k-1)*mft)
        END DO
        write(39,931) (x(j+(k-1)*nft),j=1,nft)
        write(39,931) (z(j+(k-1)*mft),j=1,mft)
       END DO
931    format(10f12.4)

  2    CALL distribution(x,z,f)
       IF((nc.eq.1).and.(nrecord1.le.1000)) toler=1.0d-06*f
       IF((nc.eq.1).and.(nrecord1.gt.1000)) toler=1.0d-05*f
       call frank
       write(43,1001) (weight(i),i=1,mft)
1001   format('Weights of features:', 20f9.5)       
       if(nc.le.1) then
        CALL clusterfunc(x,z,f)
        coef=balance*f01/f02
        print *,coef
       end if
       
       CALL rescaling(x,z,x4,z4)
       call checkfinal(x4,z4,ff,rmsetrain,rmaetrain,r2train,cortrain
     1 ,rmselarg,rmaelarg,r2larg,corlarg,rmsewei,rmaewei,r2wei,corwei
     2 ,rmseknnwei,rmaeknnwei,r2knnwei,corknnwei,rmseknn,rmaeknn,r2knn
     3 ,corknn,rmsedist,rmaedist,r2dist,cordist,rmseclust,rmaeclust
     4 ,r2clust,corclust)     

       fc(nc)=ff
       if(nc.ge.2) then
        dif1=abs(fc(nc-1)-fc(nc))/fc(1)
        if(dif1.le.1.0d-05) go to 210
       end if 

       ftr(nc,nfold)=ff
       rmsetr(nc,nfold)=rmsetrain
       rmsetest1(nc,nfold)=rmselarg
       rmsetest2(nc,nfold)=rmsewei
       rmsetest3(nc,nfold)=rmseknnwei
       rmsetest4(nc,nfold)=rmseknn
       rmsetest5(nc,nfold)=rmsedist
       rmsetest6(nc,nfold)=rmseclust

       rmaetr(nc,nfold)=rmaetrain
       rmaetest1(nc,nfold)=rmaelarg
       rmaetest2(nc,nfold)=rmaewei
       rmaetest3(nc,nfold)=rmaeknnwei
       rmaetest4(nc,nfold)=rmaeknn
       rmaetest5(nc,nfold)=rmaedist
       rmaetest6(nc,nfold)=rmaeclust

       cetr(nc,nfold)=r2train
       cetest1(nc,nfold)=r2larg
       cetest2(nc,nfold)=r2wei
       cetest3(nc,nfold)=r2knnwei
       cetest4(nc,nfold)=r2knn
       cetest5(nc,nfold)=r2dist
       cetest6(nc,nfold)=r2clust

       cortr(nc,nfold)=cortrain
       cortest1(nc,nfold)=corlarg
       cortest2(nc,nfold)=corwei
       cortest3(nc,nfold)=corknnwei
       cortest4(nc,nfold)=corknn
       cortest5(nc,nfold)=cordist
       cortest6(nc,nfold)=corclust

       WRITE(40,*)
       WRITE(40,573) nc
 573   format('No of hyperplanes:',i4)
       WRITE(40,*)
       WRITE(41,*)
       WRITE(41,*)
       DO k=1,nc
        WRITE(41,722) (x4(j+(k-1)*nft),j=1,nft)
       END DO
722    FORMAT(20f12.5)
       WRITE(40,543) ff
 543   format('Fit function value:',f24.8)
       WRITE(40,*)
       WRITE(40,541) ncount
 541   FORMAT('The number of linear regression problems solved:',i10)
       tlin1=tlin/dble(nrecord1)
       WRITE(40,*)
       WRITE(40,*)
       WRITE(40,511) tlin
 511   FORMAT('The number of linear function evaluations:',f20.0)
       WRITE(40,*)
       WRITE(40,519) tlin1
 519   FORMAT('Average number of linear function evaluations:',f20.0)
       WRITE(40,*)
       CALL cpu_time(time4)
       timef=time4-time1
       timeopt=time4-time2
       WRITE(42,610) nc,ff,ncount,tlin,tlin1,timef,timeopt
 610   FORMAT(I5,f26.6,I12,f28.0,f24.0,2f13.4)
       WRITE(40,*)
       WRITE(40,141) timef
       WRITE(40,*)
       WRITE(40,*)
       WRITE(40,*)
       IF(timef.gt.tlimit) go to 210
      END DO
141   FORMAT('CPU time:',f12.3)
      GO TO 211
 210  CONTINUE
      WRITE(43,803) 
 803  FORMAT('Fit function values on training set:')
      WRITE(43,*)
      DO i=1,nc-1
       WRITE(43,332) i
       WRITE(43,331) (ftr(i,j),j=1,nfoldmax)
      END DO
 331  FORMAT(5f28.6)
 332  FORMAT('#Clusters:',i6)
      WRITE(43,*) 
      WRITE(43,804) 
 804  FORMAT('VALUES OF PERFORMANCE MEASURES ON TRAINING SET:')       
      WRITE(43,*)      
      DO i=1,nc-1
       WRITE(43,*)
       WRITE(43,352) i
       WRITE(43,*)
       WRITE(43,351) (rmsetr(i,j),j=1,nfoldmax)
       WRITE(43,354) (rmaetr(i,j),j=1,nfoldmax)
       WRITE(43,353) (cetr(i,j),j=1,nfoldmax)
       WRITE(43,355) (cortr(i,j),j=1,nfoldmax)
      END DO
 351  FORMAT('RMSE-tr:  ',5f28.6)
 354  FORMAT('MAE-tr:   ',5f28.6)
 353  FORMAT('CD-tr:    ',5f28.6)
 355  FORMAT('PIRSON-tr:',5f28.6)
 352  FORMAT('#Clusters:',i6) 

      IF(itest.ge.2) THEN
       WRITE(43,*)
       WRITE(43,809)
 809   FORMAT('__________________________________________________')      
       WRITE(43,*)       
       WRITE(43,805) 
 805   FORMAT('VALUES OF PERFORMANCE MEASURES ON TEST SET:')         
       WRITE(43,*)      
       DO i=1,nc-1
        WRITE(43,*)
        WRITE(43,362) i
        WRITE(43,*)  
        write(43,806) 
 806    FORMAT('The first method, largest cluster:')         
        write(43,*)  
        WRITE(43,361) (rmsetest1(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest1(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest1(i,j),j=1,nfoldmax)
        WRITE(43,365) (cortest1(i,j),j=1,nfoldmax)
        write(43,*)        
        write(43,807) 
 807    FORMAT('The second method, using weights:')         
        write(43,*)  
        WRITE(43,361) (rmsetest2(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest2(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest2(i,j),j=1,nfoldmax)
        WRITE(43,365) (cortest2(i,j),j=1,nfoldmax)
        write(43,*)
        write(43,808) 
 808    FORMAT('The third method, k-NN with weights:')         
        write(43,*)  
        WRITE(43,361) (rmsetest3(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest3(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest3(i,j),j=1,nfoldmax)
        WRITE(43,365) (cortest3(i,j),j=1,nfoldmax)
        write(43,*)
        write(43,810) 
 810    FORMAT('The fourth method, k-NN one cluster:')         
        write(43,*)  
        WRITE(43,361) (rmsetest4(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest4(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest4(i,j),j=1,nfoldmax)
        WRITE(43,365) (cortest4(i,j),j=1,nfoldmax)
        write(43,*)
        write(43,811) 
 811    FORMAT('The fifth method, distances:')         
        write(43,*)  
        WRITE(43,361) (rmsetest5(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest5(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest5(i,j),j=1,nfoldmax)
        WRITE(43,365) (cortest5(i,j),j=1,nfoldmax)
        write(43,*)
        write(43,813) 
 813    FORMAT('The sixth method: cluster centers:')         
        write(43,*)  
        WRITE(43,361) (rmsetest6(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest6(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest6(i,j),j=1,nfoldmax)
        WRITE(43,365) (cortest6(i,j),j=1,nfoldmax)
       end do
 361   FORMAT('RMSE-test:  ',5f28.6)
 364   FORMAT('MAE-test:   ',5f28.6)
 363   FORMAT('CD-test:    ',5f28.6)
 365   FORMAT('PIRSON-test:',5f28.6)
 362   FORMAT('#Clusters:',i6)      
      END IF
      CLOSE(39)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(78)
      STOP 
      END
      
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Subroutine functions start here 
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>      

c======================================================================
c scaling whole data including output
c======================================================================
      SUBROUTINE scaling
      PARAMETER(maxnft=100, maxrec=200000)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      DOUBLE PRECISION a4(maxrec,maxnft),clin(maxnft),dlin(maxnft)
      COMMON /c24/a4,/crecord1/nrecord1,/cnft/nft,/ctransform/clin,dlin
c======================================================================
      DO i=1,nft
       dm=0.0d+00
       DO j=1,nrecord1
        dm=dm+a4(j,i)
       END DO
       dm=dm/DBLE(nrecord1)
       var=0.0d+00
       DO j=1,nrecord1
         var=var+(a4(j,i)-dm)**2
       END DO
       var=var/DBLE(nrecord1-1)
       var=dsqrt(var)
       IF(var.ge.1.0d-04) THEN
            clin(i)=1.0d+00/var
            dlin(i)=-dm/var
            DO j=1,nrecord1
              a4(j,i)=clin(i)*a4(j,i)+dlin(i)
            END DO
       END IF
       IF(var.lt.1.0d-04) THEN
         IF(dabs(dm).le.1.0d+00) THEN
          clin(i)=1.0d+00
          dlin(i)=0.0d+00
         END IF
         IF(dabs(dm).gt.1.0d+00) THEN
          clin(i)=1.0d+00/dm
          dlin(i)=0.0d+00
         END IF
         DO j=1,nrecord1
          a4(j,i)=clin(i)*a4(j,i)+dlin(i)
         END DO
       END IF
      END DO
cc==================================================
      RETURN
      END

c=================================================================
c  rescaling 
c==============================================================
      SUBROUTINE rescaling(x,z,x4,z4)
      PARAMETER(maxvar=1000, maxnft=100, maxrec=200000, maxclust=50)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      DOUBLE PRECISION clin(maxnft),dlin(maxnft),x(maxvar),x4(maxvar)
     1 ,a(maxrec,maxnft),z(maxvar),z4(maxvar)
      INTEGER nob(maxclust,maxrec)
     2 ,nel(maxclust)
      COMMON /cnft/nft,/cmft/mft,/ctransform/clin,dlin,/cnc/nc,/c22/a
     1 ,/crecord1/nrecord1,/cnk/nel,nob
c======================================================================
      rabs=dabs(clin(nft))
      IF(rabs.gt.1.0d-08) THEN
       DO i=1,nc
        d0=0.0d+00
        DO j=1,mft
         x4(j+(i-1)*nft)=x(j+(i-1)*nft)*clin(j)/rabs
         d0=d0+x(j+(i-1)*nft)*dlin(j)
        END DO
        x4(i*nft)=(x(i*nft)+d0-dlin(nft))/rabs
       END DO

       DO i=1,nc
        DO j=1,mft
         z4(j+(i-1)*mft)=(z(j+(i-1)*mft)-dlin(j))/clin(j)
        END DO
       END DO
      END IF   
      RETURN
      END

c======================================================================
c  distribution over clusters: here we calculate s_{k-1}^i=dminim(i)
c======================================================================
      SUBROUTINE distribution(x,z,f)
      PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      allocatable l3(:)
      DOUBLE PRECISION x(maxvar),dminim(maxrec),a4(maxrec,maxnft)
     1 ,z(maxvar),dv(maxnft),vki(maxnft),vbar(maxnft),weight(maxnft)
      INTEGER nob(maxclust,maxrec),nel(maxclust),l1(maxrec),l2(maxrec)
     1 ,nob1(maxclust,maxrec),nel1(maxclust),nob2(maxrec)
      COMMON /cnft/nft,/cmft/mft,/crecord1/nrecord1,/cnk/nel,nob,/cnc/nc
     1 ,/c24/a4,/cl1/l1,/ctlin/tlin,/cdminim/dminim,/ccoef/coef,/cl2/l2
     2 ,/ccbar/cbar,/cvbar/vbar,/cweight/weight,/cnk1/nel1,nob1
     3 ,/cnel1/nel2,nob2
      allocate(l3(maxrec))
      
      DO i=1,nc
       nel(i)=0
       nel1(i)=0
      END DO

      do j=1,mft
       vbar(j)=0.0d+00
      end do

      f=0.0d+00
      cbar=0.0d+00
      DO k=1,nrecord1
       d1=1.0d+30
       d13=1.0d+30
       d14=1.0d+30
       DO i=1,nc
        d2=x(i*nft)
        DO j=1,mft
         d2=d2+a4(k,j)*x(j+(i-1)*nft)    
        END DO
        tlin=tlin+1.0d+00
        d3=(d2-a4(k,nft))**2                                            !d3=c_{ki}
        d4=0.0d+00
        do j=1,mft
         dv(j)=dabs(z(j+(i-1)*mft)-a4(k,j))
         d4=d4+weight(j)*dv(j)
        end do
        d5=d3+coef*d4
        IF(d5.lt.d1) THEN
         d1=d5
         i1=i
         cki=d3
         do j=1,mft
          vki(j)=dv(j)
         end do
        END IF
        IF(d3.lt.d13) THEN
         d13=d3
         i3=i
        END IF
        IF(d4.lt.d14) THEN
         d14=d4
         i2=i
        END IF
       END DO
       f=f+d1
       nel(i1)=nel(i1)+1
       nob(i1,nel(i1))=k
       dminim(k)=d1
       l1(k)=i1
       l2(k)=i2
       l3(k)=i3
       cbar=cbar+cki
       do j=1,mft
        vbar(j)=vbar(j)+coef*vki(j)
       end do
      END DO
      do i=1,nrecord1
       i1=l2(i)
       nel1(i1)=nel1(i1)+1
       nob1(i1,nel1(i1))=i
      end do
      vbar1=0.0d+00
      do i=1,mft
       vbar1=dmax1(vbar1,dabs(vbar(i)))
      end do
      do i=1,mft
       vbar(i)=vbar(i)/vbar1
      end do
      cbar=cbar/vbar1
c=====================================================================
      RETURN
      END

c===================================================================
c  final distribution over clusters
c===================================================================
      subroutine checkfinal(x,z,fval,rmsetrain,rmaetrain,r2train
     1 ,cortrain,rmselarg,rmaelarg,r2larg,corlarg,rmsewei,rmaewei
     2 ,r2wei,corwei,rmseknnwei,rmaeknnwei,r2knnwei,corknnwei,rmseknn
     3 ,rmaeknn,r2knn,corknn,rmsedist,rmaedist,r2dist,cordist,rmseclust
     4 ,rmaeclust,r2clust,corclust)
      PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100)
      implicit double precision (a-h,o-z)
      allocatable ebar(:),list(:)
      double precision x(maxvar),a(maxrec,maxnft),w(maxclust)
     1 ,xw(maxnft),xc(maxvar),dist(maxclust),weight(maxnft),z(maxvar)
     2 ,bb(maxclust)
      integer nlist(maxrec),listtest(maxrec,50),nel1(maxclust)
     1 ,nob(maxclust,maxrec),nel(maxclust),nob2(maxclust,maxrec)
     2 ,nel2(maxclust)
      common /c22/a,/cnft/nft,/cnlist/nl,nlist,/cnrecord2/nrecord2
     1 ,/crecord1/nrecord1,/crecord/nrecord,/cnc/nc,/cnk1/nel,nob
     2 ,/cltest/knn,listtest,/cweight/weight,/cnk/nel2,nob2
      allocate(ebar(maxrec))
      allocate(list(maxrec))
c====================================================================
      mf=nft-1
c====================================================================      
      DO i=1,nc
       i1=nel(i)
       DO i2=1,i1
        i3=nob(i,i2)
        list(i3)=i 
       END DO
      END DO

      b0=0.0d+00
      do i=1,nc
       i1=nel(i)
       bb(i)=0.0d+00
       do j=1,i1
        j1=nob(i,j)
        b0=b0+a(j1,nft)
        bb(i)=bb(i)+a(j1,nft)
       end do
       if(nel(i).gt.0) then
           bb(i)=bb(i)/dble(nel(i))
         else
           bb(i)=0.0d+00  
       end if
      end do
      b0=b0/dble(nrecord1)

      if(nrecord2.gt.0) then
       bbar=0.0d+00  
       do i=1,nl
        i1=nlist(i)
        bbar=bbar+a(i1,nft)
       end do
       bbar=bbar/dble(nrecord2)
      end if

      res=0.0d+00
      do i=1,nc
       i1=nel(i)
       do j=1,i1
        j1=nob(i,j)
        res=res+(a(j1,nft)-b0)**2
       end do
      end do

      rtest=0.0d+00
      do i=1,nl
       i1=nlist(i) 
       rtest=rtest+(a(i1,nft)-bbar)**2
      end do
c====================================================== 
c calculation of performance measures for training set
c======================================================
      fval=0.0d+00
      rmaetrain=0.0d+00
      bmean=0.0d+00
      do i=1,nrecord
       do i1=1,nl
        IF(i.eq.nlist(i1)) GO TO 1
       end do
       f1=1.0d+26
       f2=1.0d+26
       do k=1,nc
        f3=x(k*nft)
        do j=1,nft-1
         f3=f3+a(i,j)*x(j+(k-1)*nft)
        end do
        f4=(f3-a(i,nft))**2
        f5=dabs(f3-a(i,nft))
        if(f4.lt.f1) then
         f1=f4
         ebar(i)=f3
        end if
        f2=dmin1(f2,f5)
       end do
       fval=fval+f1
       rmaetrain=rmaetrain+f2
       bmean=bmean+ebar(i)
  1   end do
      rmsetrain=sqrt(fval/dble(nrecord1))
      rmaetrain=rmaetrain/dble(nrecord1)
      r2train=1.0d+00-fval/res
      bmean=bmean/dble(nrecord1)

      r3=0.0d+00
      r4=0.0d+00
      do i=1,nrecord
       do i1=1,nl
        IF(i.eq.nlist(i1)) GO TO 2 
       end do
       r3=r3+(a(i,nft)-b0)*(ebar(i)-bmean)
       r4=r4+(ebar(i)-bmean)**2
  2   end do
      cortrain=r3/sqrt(res*r4)
c=====================================================================
c Using different methods for training set
c=====================================================================
       write(43,*)
       write(43,50) nc
 50    format('The number of clusters = ',i10)       
       write(43,*)
c======================================================================
c largest cluster 
c======================================================================      
       k1=0
       do i=1,nc
        if(k1.lt.nel(i)) then
         k1=nel(i)
         kmax=i
        end if
       end do
       rmselarg1=0.0d+00
       rmaelarg1=0.0d+00
       bmean1=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 11 
        end do
        f3=x(kmax*nft)
        do j=1,nft-1
         f3=f3+a(i,j)*x(j+(kmax-1)*nft)
        end do
        f1=(f3-a(i,nft))**2
        f2=dabs(f3-a(i,nft))
        ebar(i)=f3
        rmselarg1=rmselarg1+f1
        rmaelarg1=rmaelarg1+f2
        bmean1=bmean1+ebar(i)
 11    end do
       r2larg1=1.0d+00-rmselarg1/res
       rmselarg1=sqrt(rmselarg1/dble(nrecord1))
       rmaelarg1=rmaelarg1/dble(nrecord1)
       bmean1=bmean1/dble(nrecord1)

       r3=0.0d+00
       r4=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 12 
        end do
        r3=r3+(a(i,nft)-b0)*(ebar(i)-bmean1)
        r4=r4+(ebar(i)-bmean1)**2
  12   end do
       corlarg1=r3/sqrt(res*r4)
       write(43,44) 
44     format('First method: largest cluster')
       write(43,*)
       write(43,45) rmselarg1,rmaelarg1,r2larg1,corlarg1
45     format(4f12.6)       
       write(43,*)
c======================================================================
c The use of weights - training set
c======================================================================
       do i=1,nc
        w(i)=dble(nel(i))/dble(nrecord1)
       end do

       do j=1,nft
        xw(j)=0.0d+00
        do k=1,nc
         xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
        end do
       end do

       rmsewei1=0.0d+00
       rmaewei1=0.0d+00
       bmean1=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 13 
        end do
        f3=xw(nft)
        do j=1,nft-1
         f3=f3+a(i,j)*xw(j)
        end do
        f1=(f3-a(i,nft))**2
        f2=dabs(f3-a(i,nft))
        ebar(i)=f3
        rmsewei1=rmsewei1+f1
        rmaewei1=rmaewei1+f2
        bmean1=bmean1+ebar(i)
  13   end do
       r2wei1=1.0d+00-rmsewei1/res
       rmsewei1=sqrt(rmsewei1/dble(nrecord1))
       rmaewei1=rmaewei1/dble(nrecord1)
       bmean1=bmean1/dble(nrecord1)

       r3=0.0d+00
       r4=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 14 
        end do
        r3=r3+(a(i,nft)-b0)*(ebar(i)-bmean1)
        r4=r4+(ebar(i)-bmean1)**2
  14   end do
       corwei1=r3/sqrt(res*r4)
       write(43,46) 
46     format('Second method: weights')
       write(43,*)
       write(43,45) rmsewei1,rmaewei1,r2wei1,corwei1
       write(43,*)
c======================================================================
c The use of k nearest neighbors weighting - training set
c======================================================================
       rmseknnwei1=0.0d+00
       rmaeknnwei1=0.0d+00
       bmean1=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 15 
        end do       
        do k=1,nc
         w(k)=0.0d+00
        end do
        
        do j=1,knn
         j1=listtest(i,j)
         k=list(j1)
         w(k)=w(k)+1.0d+00/dble(knn)
        end do
       
        do j=1,nft
         xw(j)=0.0d+00
         do k=1,nc
          xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
         end do
        end do
        f3=xw(nft)
        do j=1,nft-1
         f3=f3+a(i,j)*xw(j)
        end do
        f1=(f3-a(i,nft))**2
        f2=dabs(f3-a(i,nft))
        ebar(i)=f3
        rmseknnwei1=rmseknnwei1+f1
        rmaeknnwei1=rmaeknnwei1+f2
        bmean1=bmean1+ebar(i)
 15    end do
       r2knnwei1=1.0d+00-rmseknnwei1/res
       rmseknnwei1=sqrt(rmseknnwei1/dble(nrecord1))
       rmaeknnwei1=rmaeknnwei1/dble(nrecord1)
       bmean1=bmean1/dble(nrecord1)
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 16 
        end do
        r3=r3+(a(i,nft)-b0)*(ebar(i)-bmean1)
        r4=r4+(ebar(i)-bmean1)**2
  16   end do
       corknnwei1=r3/sqrt(res*r4)
       write(43,47) 
47     format('Third method: knn with weights')
       write(43,*)
       write(43,45) rmseknnwei1,rmaeknnwei1,r2knnwei1,corknnwei1
       write(43,*)
c======================================================================
c The use of k nearest neighbors - one cluster - training set 
c======================================================================
       rmseknn1=0.0d+00
       rmaeknn1=0.0d+00
       bmean1=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 17 
        end do       
        do k=1,nc
         nel1(k)=0
        end do
        do j=1,knn
         j1=listtest(i,j)
         k=list(j1)
         nel1(k)=nel1(k)+1
        end do
        k1=0
        do i2=1,nc
         if(k1.lt.nel1(i2)) then
          k1=nel1(i2)
          kmax=i2
         end if
        end do
        f3=x(kmax*nft)
        do j=1,nft-1
         f3=f3+a(i,j)*x(j+(kmax-1)*nft)
        end do
        f1=(f3-a(i,nft))**2
        f2=dabs(f3-a(i,nft))
        ebar(i)=f3
        rmseknn1=rmseknn1+f1
        rmaeknn1=rmaeknn1+f2
        bmean1=bmean1+ebar(i)
  17   end do
       r2knn1=1.0d+00-rmseknn1/res
       rmseknn1=sqrt(rmseknn1/dble(nrecord1))
       rmaeknn1=rmaeknn1/dble(nrecord1)
       bmean1=bmean1/dble(nrecord1)
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 18 
        end do
        r3=r3+(a(i,nft)-b0)*(ebar(i)-bmean1)
        r4=r4+(ebar(i)-bmean1)**2
  18   end do
       corknn1=r3/sqrt(res*r4)
       write(43,48) 
48     format('Fourth method: knn with one cluster')
       write(43,*)
       write(43,45) rmseknn1,rmaeknn1,r2knn1,corknn1
       write(43,*)       
c======================================================================
c The use of distances - training set
c======================================================================
       do k=1,nc
        do j=1,mf
         xc(j+(k-1)*mf)=z(j+(k-1)*mf)
        end do
       end do
       rmsedist1=0.0d+00
       rmaedist1=0.0d+00
       bmean1=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 19 
        end do
        do k=1,nc
         dist(k)=0.0d+00
         if(nel(k).gt.0) then
          do j=1,mf
           dist(k)=dist(k)+(a(i,j)-xc(j+(k-1)*mf))**2
          end do
          dist(k)=sqrt(dist(k))
         end if
        end do
        
        km=0
        dc=1.0d+30
        do k=1,nc
         dc1=0.0d+00
         do k1=1,mf
          dc1=dc1+(a(i,k1)-xc(k1+(k-1)*mf))**2
         end do
         dc1=sqrt(dc1)
         if(dc.gt.dc1) then
          dc=dc1
          km=k         
         end if 
        end do
        
        if(dc.gt.1.0d-02) then
         d0=0.0d+00
         do k=1,nc
          if(nel(k).gt.0) d0=d0+1.0d+00/dist(k)**2
         end do

         do k=1,nc
          if(nel(k).gt.0) then
               w(k)=1.0d+00/(dist(k)**2*d0)
            else
               w(k)=0.0d+00  
          end if   
         end do        
        end if

        if(dc.le.1.0d-02) then
         do k=1,nc
          w(k)=0.0d+00
         end do
         w(km)=1.0d+00
        end if

        do j=1,nft
         xw(j)=0.0d+00
         do k=1,nc
          xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
         end do
        end do
        f3=xw(nft)
        do j=1,nft-1
         f3=f3+a(i,j)*xw(j)
        end do
        f1=(f3-a(i,nft))**2
        f2=dabs(f3-a(i,nft))
        ebar(i)=f3
        rmsedist1=rmsedist1+f1
        rmaedist1=rmaedist1+f2
        bmean1=bmean1+ebar(i)
  19   end do
       r2dist1=1.0d+00-rmsedist1/res
       rmsedist1=sqrt(rmsedist1/dble(nrecord1))
       rmaedist1=rmaedist1/dble(nrecord1)
       bmean1=bmean1/dble(nrecord1)

       r3=0.0d+00
       r4=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 20 
        end do
        r3=r3+(a(i,nft)-b0)*(ebar(i)-bmean1)
        r4=r4+(ebar(i)-bmean1)**2
  20   end do
       cordist1=r3/sqrt(res*r4)
       write(43,49) 
49     format('Fifth method: distance')
       write(43,*)
       write(43,45) rmsedist1,rmaedist1,r2dist1,cordist1
       write(43,*)       
c======================================================================
c Sixth method: cluster centers
c======================================================================
       rmseclust1=0.0d+00
       rmaeclust1=0.0d+00
       bmean1=0.0d+00
       DO i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 21
        end do       
        d1=1.0d+30
        do j=1,nc
         d2=0.0d+00
         do k=1,mf
          d2=d2+weight(k)*(a(i,k)-z(k+(j-1)*mf))**2
         end do
         if(d2.lt.d1) then
          d1=d2
          i2=j
         end if
        end do
        f3=x(i2*nft)
        do j=1,nft-1
         f3=f3+a(i,j)*x(j+(i2-1)*nft)
        end do
        f1=(f3-a(i,nft))**2
        f2=dabs(f3-a(i,nft))
        ebar(i)=f3
        rmseclust1=rmseclust1+f1
        rmaeclust1=rmaeclust1+f2
        bmean1=bmean1+ebar(i)
 21    END DO
       r2clust1=1.0d+00-rmseclust1/res
       rmseclust1=sqrt(rmseclust1/dble(nrecord1))
       rmaeclust1=rmaeclust1/dble(nrecord1)
       bmean1=bmean1/dble(nrecord1)
       
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nrecord
        do i1=1,nl
         IF(i.eq.nlist(i1)) GO TO 22
        end do             
        r3=r3+(a(i,nft)-b0)*(ebar(i)-bmean1)
        r4=r4+(ebar(i)-bmean)**2
 22    end do
       corclust1=r3/sqrt(res*r4)       
       write(43,51) 
51     format('Sixth method: cluster centers')
       write(43,*)
       write(43,45) rmseclust1,rmaeclust1,r2clust1,corclust1
       write(43,*)       
c======================================================================
c TEST SET
c======================================================================
c======================================================================
c Test set - the use of the largest cluster 
c======================================================================
      if(nrecord2.gt.0) then
       k1=0
       do i=1,nc
        if(k1.lt.nel(i)) then
         k1=nel(i)
         kmax=i
        end if
       end do

       rmselarg=0.0d+00
       rmaelarg=0.0d+00
       bmean=0.0d+00
       do i=1,nl
        i1=nlist(i)
        f3=x(kmax*nft)
        do j=1,nft-1
         f3=f3+a(i1,j)*x(j+(kmax-1)*nft)
        end do
        f1=(f3-a(i1,nft))**2
        f2=dabs(f3-a(i1,nft))
        ebar(i1)=f3
        rmselarg=rmselarg+f1
        rmaelarg=rmaelarg+f2
        bmean=bmean+ebar(i1)
       end do
       r2larg=1.0d+00-rmselarg/rtest
       rmselarg=sqrt(rmselarg/dble(nrecord2))
       rmaelarg=rmaelarg/dble(nrecord2)
       bmean=bmean/dble(nrecord2)

       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)       
        r3=r3+(a(i1,nft)-bbar)*(ebar(i1)-bmean)
        r4=r4+(ebar(i1)-bmean)**2
       end do
       corlarg=r3/sqrt(rtest*r4)
      end if 
c======================================================================
c Test set - the use of weights 
c======================================================================
      if(nrecord2.gt.0) then
       do i=1,nc
        w(i)=dble(nel(i))/dble(nrecord1)
       end do

       do j=1,nft
        xw(j)=0.0d+00
        do k=1,nc
         xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
        end do
       end do

       rmsewei=0.0d+00
       rmaewei=0.0d+00
       bmean=0.0d+00
       do i=1,nl
        i1=nlist(i)
        f3=xw(nft)
        do j=1,nft-1
         f3=f3+a(i1,j)*xw(j)
        end do
        f1=(f3-a(i1,nft))**2
        f2=dabs(f3-a(i1,nft))
        ebar(i1)=f3
        rmsewei=rmsewei+f1
        rmaewei=rmaewei+f2
        bmean=bmean+ebar(i1)
       end do
       r2wei=1.0d+00-rmsewei/rtest
       rmsewei=sqrt(rmsewei/dble(nrecord2))
       rmaewei=rmaewei/dble(nrecord2)
       bmean=bmean/dble(nrecord2)

       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)       
        r3=r3+(a(i1,nft)-bbar)*(ebar(i1)-bmean)
        r4=r4+(ebar(i1)-bmean)**2
       end do
       corwei=r3/sqrt(rtest*r4)
      end if 
c======================================================================
c Test set - the use of k nearest neighbors - weights 
c======================================================================
      if(nrecord2.gt.0) then
       rmseknnwei=0.0d+00
       rmaeknnwei=0.0d+00
       bmean=0.0d+00
       do i=1,nl
        i1=nlist(i)
        do k=1,nc
         w(k)=0.0d+00
        end do
        
        do j=1,knn
         j1=listtest(i,j)
         k=list(j1)
         w(k)=w(k)+1.0d+00/dble(knn)
        end do
          
        do j=1,nft
         xw(j)=0.0d+00
         do k=1,nc
          xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
         end do
        end do
        f3=xw(nft)
        do j=1,nft-1
         f3=f3+a(i1,j)*xw(j)
        end do
        f1=(f3-a(i1,nft))**2
        f2=dabs(f3-a(i1,nft))
        ebar(i1)=f3
        rmseknnwei=rmseknnwei+f1
        rmaeknnwei=rmaeknnwei+f2
        bmean=bmean+ebar(i1)
       end do
       r2knnwei=1.0d+00-rmseknnwei/rtest
       rmseknnwei=sqrt(rmseknnwei/dble(nrecord2))
       rmaeknnwei=rmaeknnwei/dble(nrecord2)
       bmean=bmean/dble(nrecord2)

       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)       
        r3=r3+(a(i1,nft)-bbar)*(ebar(i1)-bmean)
        r4=r4+(ebar(i1)-bmean)**2
       end do
       corknnwei=r3/sqrt(rtest*r4)
      end if
c======================================================================
c Test set - the use of k nearest neighbors - one cluster 
c======================================================================
      if(nrecord2.gt.0) then
       rmseknn=0.0d+00
       rmaeknn=0.0d+00
       bmean=0.0d+00
       do i=1,nl
        i1=nlist(i)
        do k=1,nc
         nel1(k)=0
        end do
        do j=1,knn
         j1=listtest(i,j)
         k=list(j1)
         nel1(k)=nel1(k)+1
        end do
        k1=0
        do i2=1,nc
         if(k1.lt.nel1(i2)) then
          k1=nel1(i2)
          kmax=i2
         end if
        end do
        f3=x(kmax*nft)
        do j=1,nft-1
         f3=f3+a(i1,j)*x(j+(kmax-1)*nft)
        end do
        f1=(f3-a(i1,nft))**2
        f2=dabs(f3-a(i1,nft))
        ebar(i1)=f3
        rmseknn=rmseknn+f1
        rmaeknn=rmaeknn+f2
        bmean=bmean+ebar(i1)
       end do
       r2knn=1.0d+00-rmseknn/rtest
       rmseknn=sqrt(rmseknn/dble(nrecord2))
       rmaeknn=rmaeknn/dble(nrecord2)
       bmean=bmean/dble(nrecord2)
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)       
        r3=r3+(a(i1,nft)-bbar)*(ebar(i1)-bmean)
        r4=r4+(ebar(i1)-bmean)**2
       end do
       corknn=r3/sqrt(rtest*r4)
      end if
c======================================================================
c Test set - the use of distances 
c======================================================================
      if(nrecord2.gt.0) then
       do k=1,nc
        do j=1,mf
         xc(j+(k-1)*mf)=0.0d+00
        end do
       end do
       do k=1,nc
         do j=1,mf
          xc(j+(k-1)*mf)=z(j+(k-1)*mf)
         end do
       end do      

       rmsedist=0.0d+00
       rmaedist=0.0d+00
       bmean=0.0d+00
       do i=1,nl
        i1=nlist(i)
        do k=1,nc
         dist(k)=0.0d+00
         if(nel(k).gt.0) then
          do j=1,mf
           dist(k)=dist(k)+(a(i1,j)-xc(j+(k-1)*mf))**2
          end do
          dist(k)=sqrt(dist(k))
         end if
        end do
        
        km=0
        dc=1.0d+30
        do k=1,nc
         dc1=0.0d+00
         do k1=1,mf
          dc1=dc1+(a(i1,k1)-xc(k1+(k-1)*mf))**2
         end do
         dc1=sqrt(dc1)
         if(dc.gt.dc1) then
          dc=dc1
          km=k         
         end if 
        end do

        if(dc.gt.1.0d-02) then
         d0=0.0d+00
         do k=1,nc
          if(nel(k).gt.0) d0=d0+1.0d+00/dist(k)**2
         end do

         do k=1,nc
          if(nel(k).gt.0) then
               w(k)=1.0d+00/(dist(k)**2*d0)
            else
               w(k)=0.0d+00  
          end if   
         end do        
        end if

        if(dc.le.1.0d-02) then
         do k=1,nc
          w(k)=0.0d+00
         end do
         w(km)=1.0d+00
        end if

        do j=1,nft
         xw(j)=0.0d+00
         do k=1,nc
          xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
         end do
        end do
        f3=xw(nft)
        do j=1,nft-1
         f3=f3+a(i1,j)*xw(j)
        end do
        f1=(f3-a(i1,nft))**2
        f2=dabs(f3-a(i1,nft))
        ebar(i1)=f3
        rmsedist=rmsedist+f1
        rmaedist=rmaedist+f2
        bmean=bmean+ebar(i1)
       end do
       r2dist=1.0d+00-rmsedist/rtest
       rmsedist=sqrt(rmsedist/dble(nrecord2))
       rmaedist=rmaedist/dble(nrecord2)
       bmean=bmean/dble(nrecord2)

       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)       
        r3=r3+(a(i1,nft)-bbar)*(ebar(i1)-bmean)
        r4=r4+(ebar(i1)-bmean)**2
       end do
       cordist=r3/sqrt(rtest*r4)
      end if
c=====================================================================
c Test set - use of cluster centers
c=====================================================================
      if(nrecord2.gt.0) then
       rmseclust=0.0d+00
       rmaeclust=0.0d+00
       bmean=0.0d+00
       DO i=1,nl
        i1=nlist(i)
        d1=1.0d+30
        do j=1,nc
         d2=0.0d+00
         do k=1,mf
          d2=d2+weight(k)*(a(i1,k)-z(k+(j-1)*mf))**2
         end do
         if(d2.lt.d1) then
          d1=d2
          i2=j
         end if
        end do
        f3=x(i2*nft)
        do j=1,nft-1
         f3=f3+a(i1,j)*x(j+(i2-1)*nft)
        end do
        f1=(f3-a(i1,nft))**2
        f2=dabs(f3-a(i1,nft))
        ebar(i1)=f3
        rmseclust=rmseclust+f1
        rmaeclust=rmaeclust+f2
        bmean=bmean+ebar(i1)
       END DO
       r2clust=1.0d+00-rmseclust/rtest
       rmseclust=sqrt(rmseclust/dble(nrecord2))
       rmaeclust=rmaeclust/dble(nrecord2)
       bmean=bmean/dble(nrecord2)
       
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)       
        r3=r3+(a(i1,nft)-bbar)*(ebar(i1)-bmean)
        r4=r4+(ebar(i1)-bmean)**2
       end do
       corclust=r3/sqrt(rtest*r4)       
      end if
c=====================================================================
c Test set - use of the mean output
c=====================================================================
      if(nrecord2.gt.0) then
       rmseout=0.0d+00
       rmaeout=0.0d+00
       bmean=0.0d+00
       DO i=1,nl
        i1=nlist(i)
        d1=1.0d+30
        do j=1,nc
         d2=x(j*nft)
         do k=1,nft-1
          d2=d2+a(i1,k)*x(k+(j-1)*nft)
         end do
         d3=abs(bb(j)-d2)
         if(d3.lt.d1) then
          d1=d3
          f3=d2
         end if
        end do
        f1=(f3-a(i1,nft))**2
        f2=dabs(f3-a(i1,nft))
        ebar(i1)=f3
        rmseout=rmseout+f1
        rmaeout=rmaeout+f2
        bmean=bmean+ebar(i1)
       END DO
       r2out=1.0d+00-rmseout/rtest
       rmseout=sqrt(rmseout/dble(nrecord2))
       rmaeout=rmaeout/dble(nrecord2)
       bmean=bmean/dble(nrecord2)
       
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)       
        r3=r3+(a(i1,nft)-bbar)*(ebar(i1)-bmean)
        r4=r4+(ebar(i1)-bmean)**2
       end do
       corout=r3/sqrt(rtest*r4)
       write(39,301) rmseout,rmaeout,r2out,corout
      end if
301   format(4f12.5)
c====================================================================
      return
      end

c====================================================================
c Step1 initialize clusters  - calculating set of initial points 
c x - is given
c xc - initial values for next linear function
c zc - initial points for the next cluster center
c====================================================================
      SUBROUTINE step1(x,xc,zc,jpoints)
      PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100
     1 ,maxinit=10000)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      ALLOCATABLE d1(:,:),d5(:)
      DOUBLE PRECISION x(maxvar),a4(maxrec,maxnft),xc(maxinit,maxnft)
     1 ,x1(maxnft),dminim(maxrec),zc(maxinit,maxnft),weight(maxnft)
      INTEGER nel(maxclust),nob(maxclust,maxrec),l1(maxrec)
      COMMON /c24/a4,/cnft/nft,/cmft/mft,/crecord1/nrecord1,/cnk/nel,nob
     1 ,/cgamma/gamma1,/cdminim/dminim,/cnc/nc,/ctoler/toler,/ctlin/tlin
     2 ,/cl1/l1,/ccoef/coef,/cweight/weight
      ALLOCATE(d1(maxclust,maxrec))
      ALLOCATE(d5(maxrec))

      DO i=1,nc-1
       DO k=1,nrecord1
        d0=0.0d+00
        DO j=1,mft
          d0=d0+x(j+(i-1)*nft)*a4(k,j)
        END DO
        tlin=tlin+1.0d+00
        d1(i,k)=d0
       END DO
      END DO

      d6=0.0d+00
      DO i=1,nrecord1
       i1=l1(i)
       DO j=1,mft
        x1(j)=x(j+(i1-1)*nft)
       END DO
       x1(nft)=a4(i,nft)-d1(i1,i)
       d5(i)=0.0d+00
       DO k=1,nrecord1
        d2=d1(i1,k)+x1(nft)
        d3=dabs(d2-a4(k,nft))
        d8=0.0d+00 
        DO j=1,mft
         d8=d8+weight(j)*dabs(a4(i,j)-a4(k,j))
        END DO
        d3=d3+coef*d8
        d4=dmin1(0.0d+00,d3-dminim(k))
        d5(i)=d5(i)+d4                                                  ! decrease of the CLR function if i-th point is taking as one cluster
       END DO
       d6=dmin1(d6,d5(i))
      END DO
      jpoints=0
      d7=gamma1*d6
      DO k=1,nrecord1
       i1=l1(k)
       IF(d5(k).le.d7) THEN
        jpoints=jpoints+1
        DO j=1,mft
         xc(jpoints,j)=x(j+(i1-1)*nft)
         zc(jpoints,j)=a4(k,j)
        END DO
        xc(jpoints,nft)=a4(k,nft)-d1(i1,k)
       END IF
      END DO
      CALL cleaning(jpoints,xc,zc)
      RETURN
      END

c=================================================================
c remove extra initial poinst using special procedure 
c=================================================================

      SUBROUTINE cleaning(jpoints,xc,zc)
      PARAMETER(maxnft=100, maxinit=10000)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      ALLOCATABLE z2(:),z3(:),z4(:),z5(:)
      DOUBLE PRECISION xc(maxinit,maxnft),zc(maxinit,maxnft)
      COMMON /cnft/nft,/cmft/mft,/ctoler/toler
      ALLOCATE(z2(maxinit))
      ALLOCATE(z3(maxinit))
      ALLOCATE(z4(maxinit))
      ALLOCATE(z5(maxinit))
      
      jpoints0=0
      DO i=1,jpoints
       DO j=1,nft
        z2(j)=xc(i,j)
       END DO
       DO j=1,mft
        z3(j)=zc(i,j)
       END DO
       IF(i.eq.1) THEN
        jpoints0=1
        DO j=1,nft
         z4(j)=z2(j)
        END DO
        DO j=1,mft
         z5(j)=z3(j)
        END DO
       END IF

       IF(i.gt.1) THEN
        DO k=1,jpoints0
         d2=0.0d+00
         DO j=1,nft
          d2=d2+ABS(z2(j)-z4(j+(k-1)*nft))
         END DO
         IF(d2.le.toler) GO TO 1
         d3=0.0d+00
         DO j=1,nft
          d3=d3+ABS(z3(j)-z5(j+(k-1)*nft))
         END DO
         IF(d3.le.toler) GO TO 1
        END DO
        jpoints0=jpoints0+1
        DO j=1,nft
         z4(j+(jpoints0-1)*nft)=z2(j)
        END DO
        DO j=1,mft
         z5(j+(jpoints0-1)*mft)=z3(j)
        END DO
       END IF
   1  END DO
      DO k=1,jpoints0
       DO j=1,nft
        xc(k,j)=z4(j+(k-1)*nft)
       END DO
       DO j=1,mft
        zc(k,j)=z5(j+(k-1)*nft)
       END DO
      END DO
      jpoints=jpoints0
      END SUBROUTINE

c=====================================================
      SUBROUTINE optim(x,z,fvalue)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxnft=100)
      DOUBLE PRECISION x(maxvar),z(maxvar)
      COMMON /cnft/nft,/cns/ns,/cnc/nc,/cng/ng,/csize/m,/csize1/m1
c=======================================================
c  Input data:
c  n        - number of variables
c=======================================================
      IF(ns.eq.1) THEN
        n = 2*nft-1
        m1 = nft
      END IF
      IF(ns.eq.2) THEN
        n = nc*(2*nft-1)
        m1 = nft*nc
      END IF
      m=n
c======================================================
c ng = 0 if you use approximation to subgradient
c ng = 1 if you use exact gradients of the objective
c               and constraint functions
c======================================================
      ng=1
      CALL optimum(x,z,fvalue)
c=======================================================
      RETURN
      END
c=====================================================
!  optimization method
c=====================================================
      SUBROUTINE optimum(x,z,f2)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxdg=1000, maxit=100000, maxnft=100)
      ALLOCATABLE fvalues(:),prod(:,:),w(:,:)
      DOUBLE PRECISION x(maxvar),g(maxvar),v(maxvar),z1(maxdg)
     1 , z(maxvar), y(maxvar), y1(maxvar)
      INTEGER ij(maxdg)
      COMMON /csize/m,/cij/ij,jvertex,/cz/z1,/ckmin/kmin,/cns/ns
     1 ,/csize1/m1 
      ALLOCATE(fvalues(maxit))
      ALLOCATE(prod(maxdg,maxdg))
      ALLOCATE(w(maxdg,maxvar))
c====================================================================
      dist1=1.0d-07
      step0=-2.0d-01
      div=1.0d-01
      eps0=1.0d-07
      slinit=1.0d+00
      slmin=1.0d-05*slinit
      maxiter=5000
      niter=0
      pwt=1.0d-07
      sdif=1.0d-04
      mturn=4
      nbundle=MIN(20,2*m+1)
c====================================================================
      do i=1,m1
       y(i)=x(i)
      end do
      do i=m1+1,m
       y(i)=z(i-m1)
      end do
      sl=slinit/div
      CALL func(y,f2)
  1   sl=div*sl
      IF(sl.lt.slmin) RETURN
      DO i=1,m
       g(i)=1.0d+00/dsqrt(DBLE(m))
      END DO
      nnew=0
c================================================================
   2  niter=niter+1
      f1=f2
        IF(niter.gt.maxiter) RETURN
        nnew=nnew+1
        fvalues(niter)=f1
c---------------------------------------------------------------
        IF(nnew.gt.mturn) THEN
         mturn2=niter-mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.sdif) GO TO 1
        END IF
        IF(nnew.GE.(2*mturn)) THEN
         mturn2=niter-2*mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.(1.0d+01*sdif)) GO TO 1
        END IF
c--------------------------------------------------------------
        DO ndg=1,nbundle
            CALL dgrad(y,sl,g,v,f4,ndg,pwt)
            dotprod=0.0d+00
            DO i=1,m
             dotprod=dotprod+v(i)*v(i)
            END DO
            r=dsqrt(dotprod)
            IF(r.lt.eps0) GO TO 1
            IF(ndg.eq.1) THEN
                         rmean=r
                         kmin=1
                         rmin=r
            END IF
            IF(ndg.gt.1) THEN
                         rmin=dmin1(rmin,r)
                         IF(r.eq.rmin) kmin=ndg
                         rmean=((ndg-1)*rmean+r)/ndg
            END IF
            toler=dmax1(eps0,dist1*rmean)
            DO i=1,ndg-1
             prod(ndg,i)=0.0d+00
             DO j=1,m
              prod(ndg,i)=prod(ndg,i)+w(i,j)*v(j)
             END DO
             prod(i,ndg)=prod(ndg,i)
            END DO
            prod(ndg,ndg)=dotprod
c====================================================================
            DO i=1,m
             w(ndg,i)=v(i)
            END DO
            CALL wolfe(ndg,prod)
c================================
            DO i=1,m
             v(i)=0.0d+00
             DO j=1,jvertex
              v(i)=v(i)+w(ij(j),i)*z1(j)
             END do
            END do
c================================
            r=0.0d+00
            DO i=1,m
             r=r+v(i)*v(i)
            END DO
            r=dsqrt(r)
            IF(r.lt.toler) GO TO 1
c===========================================================
             DO i=1,m
              g(i)=-v(i)/r
              y1(i)=y(i)+sl*g(i)
             END DO
c===========================================================
             CALL func(y1,f4)
             f3=(f4-f1)/sl
             decreas=step0*r
             IF(f3.lt.decreas) THEN
                        CALL armijo(y,g,f1,f5,f4,sl,step,r)
                        f2=f5
                        DO i=1,m
                         y(i)=y(i)+step*g(i)
                        END DO
                        do i=1,m1
                         x(i)=y(i)
                        end do
                        do i=m1+1,m
                         z(i-m1)=y(i)
                        end do                        
                        sl=1.2d+00*sl
                        GO TO 2
             END IF
         END do
c=====================================================
      go to 1
      RETURN
      END

c==============================================================
c  Subroutines Wolfe and Equations solves quadratic
c  programming problem, to find
c  descent direction, Step 3, Algorithm 2.
c===============================================================

      SUBROUTINE wolfe(ndg,prod)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxdg=1000)
      COMMON /csize/m,/w01/a,/cij/ij,jvertex,/cz/z,/ckmin/kmin
      INTEGER ij(maxdg)
      DOUBLE PRECISION z(maxdg),z1(maxdg),a(maxdg,maxdg)
     1 ,prod(maxdg,maxdg)
      j9=0
      jmax=500*ndg
      jvertex=1
      ij(1)=kmin
      z(1)=1.0d+00
c=======================================
c  To calculate norm of X
c=======================================
 1    r=0.0d+00
      DO i=1,jvertex
       DO j=1,jvertex
        r=r+z(i)*z(j)*prod(ij(i),ij(j))
       END DO
      END DO
      IF(ndg.eq.1) RETURN
c========================================
c  To calculate <X,P_J> and J
c========================================
      t0=1.0d+30
      DO i=1,ndg
        t1=0.0d+00
        DO j=1,jvertex
          t1=t1+z(j)*prod(ij(j),i)
        END DO
        IF(t1.lt.t0) THEN
                     t0=t1
                     kmax=i
        END IF
      END DO
c========================================
c  First stopping criterion
c========================================
      rm=prod(kmax,kmax)
      DO j=1,jvertex
       rm=dmax1(rm,prod(ij(j),ij(j)))
      END DO
      r2=r-1.0d-12*rm
      IF(t0.gt.r2) RETURN
c========================================
c  Second stopping criterion
c========================================
      DO i=1,jvertex
       IF(kmax.eq.ij(i)) RETURN
      END DO
c========================================
c Step 1(e) from Wolfe's algorithm
c========================================
      jvertex=jvertex+1
      ij(jvertex)=kmax
      z(jvertex)=0.0d+00
c========================================
 2    DO i=1,jvertex
       DO j=1,jvertex
        a(i,j)=1.0d+00+prod(ij(i),ij(j))
       END DO
      END DO
      j9=j9+1
      IF(j9.gt.jmax) RETURN
      CALL equations(jvertex,z1)
      DO i=1,jvertex
       IF(z1(i).le.1.0d-10) go to 3
      END DO
      DO i=1,jvertex
       z(i)=z1(i)
      END DO
      go to 1
  3   teta=1.0d+00
      DO i=1,jvertex
       z5=z(i)-z1(i)
       IF(z5.gt.1.0d-10) teta=dmin1(teta,z(i)/z5)
      END DO
      kzero=0
      DO i=1,jvertex
       z(i)=(1.0d+00-teta)*z(i)+teta*z1(i)
       IF(z(i).le.1.0d-10) THEN
                          z(i)=0.0d+00
                          kzero=i
       END IF
      END DO
      j2=0
      DO i=1,jvertex
       IF(i.ne.kzero) THEN
                     j2=j2+1
                     ij(j2)=ij(i)
                     z(j2)=z(i)
       END IF
      END DO
      jvertex=j2
      go to 2
      RETURN
      END

c=====================================================
c Solving systems of equations form QP problems
c=====================================================
      SUBROUTINE equations(n,z1)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxdg=1000)
      ALLOCATABLE b(:,:)
      COMMON /w01/a
      DOUBLE PRECISION a(maxdg,maxdg),z1(maxdg)
      ALLOCATE(b(maxdg,maxdg))
      DO i=1,n
       DO j=1,n
        b(i,j)=a(i,j)
       END DO
       b(i,n+1)=1.0d+00
      END DO
      DO i=1,n
       r=b(i,i)
       DO j=i,n+1
        b(i,j)=b(i,j)/r
       END DO
       DO j=i+1,n
        DO k=i+1,n+1
         b(j,k)=b(j,k)-b(i,k)*b(j,i)
        END DO
       END DO
      END DO
      z1(n)=b(n,n+1)
      DO i=1,n-1
        k=n-i
        z1(k)=b(k,n+1)
        DO j=k+1,n
         z1(k)=z1(k)-b(k,j)*z1(j)
        END do
      END DO
      z2=0.0d+00
      DO i=1,n
       z2=z2+z1(i)
      END DO
      DO i=1,n
       z1(i)=z1(i)/z2
      END DO
      RETURN
      END

c=====================================================================
c Subroutine dgrad calculates subgradients or discrete gradients
c=====================================================================
      SUBROUTINE dgrad(y,sl,g,dg,f4,ndg,pwt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxdg=1000)
      DOUBLE PRECISION y1(maxvar),g(maxvar),y(maxvar),dg(maxvar)
     1 ,x1(maxvar), z1(maxvar) 
      COMMON /csize/m,/cng/ng,/csize1/m1
     
      DO k=1,m
        y1(k)=y(k)+sl*g(k)
      END DO

      IF(ng.eq.0) THEN
       IF(ndg.gt.1) r2=f4
       IF(ndg.eq.1) CALL func(y1,r2)
       CALL dgrad2(y1,dg,r2,pwt)
      END IF

      IF(ng.eq.1) then
       do i=1,m1
        x1(i)=y1(i)
       end do
       do i=m1+1,m
        z1(i-m1)=y1(i)
       end do
       CALL fgrad(x1,z1,dg)
      end if
      RETURN
      END

c=====================================================================
c Subroutine dgrad calculates discrete gradients
c=====================================================================
      SUBROUTINE dgrad2(y1,v,r2,pwt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxdg=1000)
      DOUBLE PRECISION y1(maxvar),v(maxvar)
      COMMON /csize/m
      t=pwt
      r4=r2
      DO k=1,m
        r3=r4
        y1(k)=y1(k)+t
        CALL func(y1,r4)
        v(k)=(r4-r3)/t
      END DO
      RETURN
      END

c===========================================================
c Line search (Armijo-type)
c===========================================================
      SUBROUTINE armijo(y,g,f1,f5,f4,sl,step,r)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000)
      COMMON /csize/m
      DOUBLE PRECISION y(maxvar),g(maxvar),y1(maxvar)
      step=sl
      f5=f4
      s1=sl
      k=0
  1   k=k+1
      IF(k.gt.20) RETURN
      s1=2.0d+00*s1
      DO i=1,m
       y1(i)=y(i)+s1*g(i)
      END DO
      CALL func(y1,f50)
      f30=f50-f1+5.0d-02*s1*r
      IF(f30.gt.0.0d+00) RETURN
      step=s1
      f5=f50
      GO TO 1
      RETURN
      END

c=====================================================
!   
c=====================================================      
      SUBROUTINE func(y,objf)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000)
      DOUBLE PRECISION y(maxvar),x(maxvar),z(maxvar)
      COMMON /cns/ns,/csize/m,/csize1/m1
c===================================================================
      do i=1,m1
       x(i)=y(i)
      end do
      do i=m1+1,m
       z(i-m1)=y(i)
      end do
      IF(ns.eq.1) CALL auxfunc(x,z,f)
      IF(ns.eq.2) CALL clusterfunc(x,z,f)
      objf=f
c===================================================================
      RETURN
      END

c==================================================================
! define weighted clustering CLR function (WC-CLR)
c==================================================================
      SUBROUTINE clusterfunc(x,z,f)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100)
      DOUBLE PRECISION x(maxvar),z(maxvar),a4(maxrec,maxnft)
     1 ,weight(maxnft)
      integer km(maxrec),kmh(maxrec)
      COMMON /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/cmft/mft
     1 ,/ctlin/tlin,/ckmin1/km,/ccoef/coef,/cf1f2/f1,f2,/cweight/weight
     2 ,/cdelta/delta,/ckmh/kmh

      f=0.0d+00
      f1=0.0d+00
      f2=0.0d+00
      DO i=1,nrecord1
       d5=1.0d+30
       DO k=1,nc
        d4=0
        DO j=1,mft
          d4=d4+weight(j)*dabs(z(j+(k-1)*mft)-a4(i,j))                  ! weight*(x-a)^2 
        END DO
        d1=x(k*nft)                                                     ! y
        DO j=1,mft
         d1=d1+a4(i,j)*x(j+(k-1)*nft)                                   ! <a,x>+y
        END DO
        tlin=tlin+1.0d+00
        d61=d1-a4(i,nft)
        d6=dabs(d61)
        d2=d6**2
        if(d6.le.delta) then
            d7=5.0d-01*d2
          else
            d7=delta*(d6-5.0d-01*delta) 
        end if
        d3=coef*d4+d7
        if(d3.lt.d5) then
          km(i)=k
          d5=d3
          d51=d4
          d52=d2
          if(d6.le.delta) then
              kmh(i)=1
            else
              kmh(i)=2
          end if
        end if
       END DO
       f1=f1+d51
       f2=f2+d52
       f=f+d5
      END DO
c=====================================================
      RETURN
      END

c=================================================================
! define aux weighted clustering CLR function (aux WC-CLR)
c=================================================================
      SUBROUTINE auxfunc(x,z,fval)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100)
      DOUBLE PRECISION x(maxvar),a4(maxrec,maxnft),dminim(maxrec)
     1 ,z(maxvar),weight(maxnft)
      integer kmaux(maxrec)
      COMMON /c24/a4,/crecord1/nrecord1,/cnft/nft,/cmft/mft,/ctlin/tlin
     1 ,/cdminim/dminim,/ckmaux/kmaux,/ccoef/coef,/cweight/weight
     
      fval=0.0d+00
      DO i=1,nrecord1
       f1=dminim(i)                                                     ! s_(k-1)^i 
       f2=x(nft)                                                        ! v intercept
       DO j=1,mft
        f2=f2+a4(i,j)*x(j)                                              ! <u,a>
       END DO
       tlin=tlin+1.0d+00
       f3=(a4(i,nft)-f2)**2                                             ! (b_i-<u,a>-v)^2    CLR 
       f4=0.0d+00
       DO j=1,mft
        f7=z(j)-a4(i,j)
        f4=f4+weight(j)*dabs(f7)                                        ! weight*(q-a)^2  Weighted clust
       END DO
       f5=f3+coef*f4                                                    ! clr + weighted clust
       f6=dmin1(f1,f5)
       if(f6.eq.f5) then
          kmaux(i)=1
        else
          kmaux(i)=0
       end if
       fval=fval+f6
      END DO
      RETURN
      END

c===================================================== 
      SUBROUTINE fgrad(x,z,subg)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100)
       DOUBLE PRECISION x(maxvar), a4(maxrec,maxnft), subg(maxvar)
     1 ,z(maxvar),weight(maxnft)
       integer km(maxrec),kmaux(maxrec),kmh(maxrec)
       COMMON /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/cmft/mft
     1 ,/cns/ns,/ctlin/tlin,/csize/m,/ckmin1/km,/ckmaux/kmaux
     2 ,/ccoef/coef,/cweight/weight,/cdelta/delta,/ckmh/kmh
      
      DO i=1,m
       subg(i)=0.0d+00
      END DO

      IF(ns.eq.1) THEN
       call auxfunc(x,z,f)
       DO i=1,nrecord1
        if(kmaux(i).eq.1) then
         d1=x(nft)
         DO j=1,mft
          d1=d1+a4(i,j)*x(j)
         END DO
         tlin=tlin+1.0d+00
         d1=d1-a4(i,nft)
         DO j=1,mft
          subg(j)=subg(j)+2.0d+00*d1*a4(i,j)  
         END DO
         subg(nft)=subg(nft)+2.0d+00*d1
         do j=1,mft
          d2=z(j)-a4(i,j)
          if(d2.ge.0.0d+00) then
               subg(nft+j)=subg(nft+j)+weight(j)*coef
             else
               subg(nft+j)=subg(nft+j)-weight(j)*coef
          end if
         end do          
        end if
       END DO
      END IF

      IF(ns.eq.2) THEN
       call clusterfunc(x,z,f)
       DO i=1,nrecord1
        jmin=km(i)
        d4=x(jmin*nft)
        do j=1,mft
         d4=d4+x(j+(jmin-1)*nft)*a4(i,j)
        end do
        d4=d4-a4(i,nft)
c----------------------------------------------------------------------
        if(kmh(i).eq.1) then
         do j=1,mft
          subg(j+(jmin-1)*nft)=subg(j+(jmin-1)*nft)+d4*a4(i,j)
         END DO
         subg(jmin*nft)=subg(jmin*nft)+d4
        end if        
        
        if(kmh(i).eq.2) then
         if(d4.ge.0.0d+00) then
           do j=1,mft
            subg(j+(jmin-1)*nft)=subg(j+(jmin-1)*nft)+delta*a4(i,j)
           END DO
           subg(jmin*nft)=subg(jmin*nft)+delta
          else 
           do j=1,mft
            subg(j+(jmin-1)*nft)=subg(j+(jmin-1)*nft)-delta*a4(i,j)
           END DO
           subg(jmin*nft)=subg(jmin*nft)-delta
         end if
        end if 
c----------------------------------------------------------------------
        n3=(jmin-1)*mft
        n1=nc*nft+n3
        do j=1,mft
         d3=z(n3+j)-a4(i,j)
         if(d3.gt.0.0d+00) then
               subg(n1+j)=subg(n1+j)+weight(j)*coef
             else
               subg(n1+j)=subg(n1+j)-weight(j)*coef
         end if
        end do
       END DO
      END IF
      END SUBROUTINE 

c=====================================================
!  BFGS for one CLR
c=====================================================
      SUBROUTINE bfgs(u)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxnft=100)
      DOUBLE PRECISION gprev(maxvar),gcur(maxvar),s(maxvar)
     1 ,sk(maxvar),yk(maxvar),h(maxvar,maxvar),u(maxvar)
      COMMON /matrices/sk,yk,h,/csize2/m,/ceps/eps,/cnft/nft
       m=nft
       eps=1.0d-05
       k=0
       DO i=1,m
        DO j=1,m
         IF(i.eq.j) h(i,j)=1.0d+00
         IF(i.ne.j) h(i,j)=0.0d+00
        END DO
       END DO
       CALL funcb(u,f2)
   1   f1=f2
       CALL gradient(u,gcur)
       k=k+1
       IF(k.gt.1000) RETURN

       gradnorm=0.0d+00
       DO i=1,m
        gradnorm=gradnorm+gcur(i)**2
       END DO
       gradnorm=dsqrt(gradnorm)
       IF(gradnorm.le.eps) RETURN
       IF(k.gt.1) THEN
        d4=0.0d+00
        DO i=1,m
         yk(i)=gcur(i)-gprev(i)
         d4=d4+yk(i)**2
        END DO
        d4=SQRT(d4)
        IF(d4.le.eps) RETURN
       END IF
       CALL dir(s,gcur,k)

       d=0.0d+00
       DO i=1,m
        d=d+s(i)**2
       END DO
       d=SQRT(d)
       IF(d.le.eps) RETURN

       DO i=1,m
        s(i)=-s(i)/d
       END DO
       d2=0.0d+00
       DO i=1,m
        d2=d2-s(i)*gcur(i)
       END DO

       IF(d2.le.0.0d+00) RETURN
       du=0.0d+00
       dg=0.0d+00
       DO i=1,m
        du=du+u(i)**2
        dg=dg+u(i)*s(i)
       END DO

       stepmax=gradnorm
       CALL armijob(f1,f5,s,u,c,stepmax,d2)
       IF(c.le.eps) RETURN
       d3=0.0d+00
       DO i=1,m
        sk(i)=c*s(i)
        u(i)=u(i)+sk(i)
        gprev(i)=gcur(i)
        d3=d3+sk(i)**2
       END DO
       d3=SQRT(d3)
       IF(d3.le.eps) RETURN
       f2=f5
       dif1=ABS(f1-f2)/(ABS(f2)+1.0d+00)
       IF(dif1.le.eps) RETURN
       GO TO 1
      RETURN
      END

c=====================================================================
! Error function for one cluster and one linear function
c=====================================================================
      SUBROUTINE funcb(x,f)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxvar=1000, maxrec=200000, maxclust=50, maxnft=100)
      DOUBLE PRECISION a4(maxrec,maxnft),x(maxvar)
      INTEGER nob1(maxrec)
      COMMON /cnel1/nel1,nob1,/c24/a4,/cnft/nft,/cmft/mft,/ctlin/tlin
c=====================================================================
      f=0.0d+00
      DO k=1,nel1
        d1=x(nft)
        DO j=1,mft
         d1=d1+a4(nob1(k),j)*x(j)
        END DO
        tlin=tlin+1.0d+00
        d2=(d1-a4(nob1(k),nft))**2 
        f=f+d2
      END DO
c=====================================================================
      RETURN
      END      

c=======================================================================
! Gradient of the function for updating one line and one cluster center
c=======================================================================
      SUBROUTINE gradient(x,grad)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxrec=200000, maxclust=50, maxnft=100)
      DOUBLE PRECISION x(maxnft),a4(maxrec,maxnft),grad(maxnft)
      INTEGER nob1(maxrec)
      COMMON /cnel1/nel1,nob1,/c24/a4,/cnft/nft,/cmft/mft,/ctlin/tlin
c=====================================================================
      DO j=1,nft
       grad(j)=0.0d+00
      END DO
      DO k=1,nel1
       d1=x(nft)
       DO j=1,mft
        d1=d1+a4(nob1(k),j)*x(j)
       END DO
       tlin=tlin+1.0d+00
       d1=d1-a4(nob1(k),nft)
       DO j=1,mft
        grad(j)=grad(j)+2.0d+00*d1*a4(nob1(k),j)   
       END DO
       grad(nft)=grad(nft)+2.0d+00*d1
      END DO
c===============================================================
      RETURN
      END

c=====================================================
       SUBROUTINE dir(s,gcur,k)
       PARAMETER(maxvar=1000)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       DOUBLE PRECISION s(maxvar),gcur(maxvar),h(maxvar,maxvar)
     1  ,sk(maxvar),yk(maxvar)
       COMMON /matrices/sk,yk,h,/csize2/m
       IF(k.gt.1) CALL matrix1
       DO i=1,m
        s(i)=0.0d+00
        DO j=1,m
         s(i)=s(i)+h(i,j)*gcur(j)
        END DO
       END DO
      RETURN
      END
 
c=====================================================
      SUBROUTINE matrix1
      PARAMETER(maxvar=1000)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      ALLOCATABLE h1(:,:),h2(:,:),h3(:,:),h4(:,:)
      DOUBLE PRECISION yk(maxvar),sk(maxvar),h(maxvar,maxvar)
      COMMON /matrices/sk,yk,h,/csize2/m
       ALLOCATE(h1(maxvar,maxvar))
       ALLOCATE(h2(maxvar,maxvar))
       ALLOCATE(h3(maxvar,maxvar))
       ALLOCATE(h4(maxvar,maxvar))
       rho=0.0d+00
       DO i=1,m
        rho=rho+sk(i)*yk(i)
       END DO
       DO i=1,m
        DO j=1,m
         h1(i,j)=(sk(i)*yk(j))/rho
         h2(i,j)=(yk(i)*sk(j))/rho
         h3(i,j)=(sk(i)*sk(j))/rho
        END DO
       END DO
       DO i=1,m
        DO j=1,m
         IF(i.eq.j) h1(i,j)=1.0d+00-h1(i,j)
         IF(i.ne.j) h1(i,j)=-h1(i,j)
         IF(i.eq.j) h2(i,j)=1.0d+00-h2(i,j)
         IF(i.ne.j) h2(i,j)=-h2(i,j)
        END DO
       END DO
       DO i=1,m
        DO j=1,m
         h4(i,j)=0.0d+00
         DO k=1,m
          h4(i,j)=h4(i,j)+h(i,k)*h2(k,j)
         END DO
        END DO
       END DO
       DO i=1,m
        DO j=1,m
         h(i,j)=h3(i,j)
         DO k=1,m
          h(i,j)=h(i,j)+h1(i,k)*h4(k,j)
         END DO
        END DO
       END DO
      RETURN
      END

c=====================================================
       SUBROUTINE armijob(f1,f5,s,u,step,stepmax,r)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       PARAMETER(maxvar=1000)
       DOUBLE PRECISION u(maxvar),u1(maxvar),s(maxvar)
       COMMON /csize2/m,/ceps/eps
       k=0
       step=2.0d+00*stepmax
  1    step=5.0d-01*step
       k=k+1
       if(k.gt.20) RETURN
       DO i=1,m
        u1(i)=u(i)+step*s(i)
       END DO
       CALL funcb(u1,f5)
       f3=f5-f1+1.0d-02*step*r
       IF(f3.lt.0.0d+00) RETURN
       IF(step.lt.eps) RETURN
       GO TO 1
       RETURN
      END

c=============================================================
c Weight optimization: Frank-Wolfe method
c=====================================================================
      SUBROUTINE frank
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(maxnft=100)
      allocatable point(:,:)
      DOUBLE PRECISION w(maxnft),vbar(maxnft),dir(maxnft),d(maxnft)
     1 ,grad(maxnft) 
      COMMON /cmft/mft,/ccbar/cbar,/cvbar/vbar,/cweight/w,/cstab/tstab
      allocate(point(maxnft,maxnft)) 
c===================================================================
      t=tstab
      delta=-1.0d-04
      do k=1,mft
       do j=1,mft
        point(k,j)=0.0d+00
       end do
       point(k,k)=1.0d+00
      end do
      do k=1,mft
       w(k)=1.0d+00/dble(mft)
      end do
  1   continue
      do k=1,mft
       grad(k)=vbar(k)+t*w(k)
      end do
      d2=1.0d+30
      do k=1,mft
       do j=1,mft
        d(j)=point(k,j)-w(j)
       end do
       d1=0.0d+00
       do j=1,mft
        d1=d1+grad(j)*d(j)
       end do
       if(d1.le.d2) then
        d2=d1
        do j=1,mft
         dir(j)=d(j)
        end do
       end if
      end do
      if(d2.ge.delta) return
      d3=0.0d+00
      do j=1,mft
       d3=d3+dir(j)**2
      end do
      sigma=-d2/(t*d3)
      sigma=dmin1(1.0d+00,sigma)
      do j=1,mft
       w(j)=w(j)+sigma*dir(j)
      end do
      go to 1 
      RETURN
      END      
