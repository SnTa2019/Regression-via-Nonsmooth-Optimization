c=======================================================================
c  ADC-CLR - Adaptive DC clusterwise linear regression  
c=======================================================================
c     main programm
c=======================================================================
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100
     1 ,maxinit=1000000)
      implicit double precision(a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),a1(maxnft)
     1 ,dminim(maxrec),z1(maxdim),xc(maxinit),z(maxdim)
     2 ,fv2(maxrec),xbest(maxdim),x1(maxdim),x4(maxdim)
     3 ,a4(maxrec,maxnft),clin(maxnft),dlin(maxnft)
     3 ,rmsetr(maxclust,100),cetr(maxclust,100),ftr(maxclust,100)
     4 ,rmsetest2(maxclust,100),rmsetest3(maxclust,100)
     5 ,rmsetest4(maxclust,100),cetest2(maxclust,100)
     6 ,cetest3(maxclust,100),cetest4(maxclust,100),pirs1(maxclust,100)
     7 ,pirs2(maxclust,100),pirs3(maxclust,100),pirs4(maxclust,100)
     8 ,pirs5(maxclust,100),rmaetr(maxclust,100),rmaetest2(maxclust,100)
     9 ,rmaetest3(maxclust,100),rmaetest4(maxclust,100)
     9 ,cetest5(maxclust,100),rmsetest5(maxclust,100)
     9 ,rmaetest5(maxclust,100),d3(maxrec)    
      integer nob(maxclust,maxrec),nel(maxclust),nob1(maxrec)
     1 ,nlist(maxrec),listtest(100000,50)
      common /c22/a,/anclust/nclust,/cnft/nft,/crecord/nrecord,/c24/a4
     1 ,/cnk/nel,nob,/ctoler/toler,/cgamma/gamma1,/cdminim/dminim
     2 ,/cnc/nc,/cnel1/nel1,nob1,/cnlist/nl,nlist,/ctlin/tlin
     3 ,/crecord1/nrecord1,/cnrecord2/nrecord2,/ctransform/clin,dlin
     4 ,/cltest/knn,listtest
      open(40,file='hyperplanes.txt')
      open(42,file='statistics.txt')
      open(41,file='predictions.txt')
      open(43,file='trainingtest.txt')
      open(78,file='datainput.txt',status='old',form='formatted')      
c=======================================================================
      PRINT *,' '
      PRINT *,'Number of features (including output column):'
      read *,nft
      PRINT *,' '
      PRINT *,'Output column:'
      read *,noutcom
      PRINT *,' '
      PRINT *,'Number of clusters:'
      read *,nclust
c=======================================================================      
      do i=1,maxrec
       read(78,*,END=901,ERR=900) (a(i,k),k=1,nft)
       nrecord=i
      end do
  900 stop 'Error in input file'       
  901 WRITE(40,*) 'Input complete. Number of records: ',nrecord
      WRITE(40,*)
c=======================================================================
      PRINT *,' '
      PRINT *,'Choose:'
      PRINT *,'1: if you do not use test set'
      PRINT *,'2: if you use one training set and one test set'
      PRINT *,'3: if you use cross-validation'
      read *,itest
      if(itest.eq.1) then
        ntrain=nrecord
        nfoldmax=1
      end if
      if(itest.eq.2) then
        PRINT *,'Enter the number of points in training set'
        read *,ntrain
        nfoldmax=1
      end if
      if(itest.eq.3) then
        PRINT *,'Enter the number of cross validations:'
        read *,nfoldmax
      end if
      PRINT *,' '
c=======================================================================
      tlimit=7.2d+04
      call cpu_time(time1)
c=======================================================================   
      if(noutcom.lt.nft) then   
       do i=1,nrecord
        j1=0
        do j=1,nft
         if (j.ne.noutcom) then
          j1=j1+1
          a1(j1)=a(i,j)
         end if
        end do
        a1(nft)=a(i,noutcom)
        do j=1,nft
         a(i,j)=a1(j)
        end do
       end do
      end if
      
      if(nrecord.le.200) then
       gamma1=0.0d+00
       gamma2=2.0d+01
      end if

      if((nrecord.gt.200).and.(nrecord.le.1000)) then
       gamma1=4.5d-01
       gamma2=1.5d+00
      end if
      
      if((nrecord.gt.1000).and.(nrecord.le.4000)) then
       gamma1=8.5d-01
       gamma2=1.15d+00
      end if

      if((nrecord.gt.4000).and.(nrecord.le.15000)) then
       gamma1=9.95d-01
       gamma2=1.05d+00
      end if
      if((nrecord.gt.15000).and.(nrecord.le.50000)) then
       gamma1=9.99d-01
       gamma2=1.02d+00
      end if
      if(nrecord.gt.50000) then
       gamma1=9.999d-01
       gamma2=1.005d+00
      end if
c=======================================================================
       WRITE(42,701)
 701   FORMAT('#Clust','             fvalue','        #Regprob_eval',
     1 '           #Lin_func_val','     #Lin_eval_per_point',
     2 '          CPU')
       WRITE(42,*)
c=======================================================================
c Division to training and test sets (not randomly)
c=======================================================================
      nfold=0
  211 nfold=nfold+1
      IF(nfold.gt.nfoldmax) GO TO 210
      if(itest.le.2) nrecord2=nrecord-ntrain
      if(itest.eq.3) nrecord2=nrecord/nfoldmax
      if(itest.eq.3) then
       nrecord3=(nfold-1)*nrecord2+1
       nrecord4=nfold*nrecord2
       nl=0
       do j=nrecord3,nrecord4
        nl=nl+1
        nlist(j-nrecord3+1)=j 
       end do
       nl1=0
       do i=1,nrecord
        do j=1,nl
         IF(i.eq.nlist(j)) GO TO 203
        end do
        nl1=nl1+1
        do k=1,nft
         a4(nl1,k)=a(i,k)
        end do
 203   end do
      end if
      if(itest.eq.1) then
       nl=0
       do i=1,nrecord
        do k=1,nft
         a4(i,k)=a(i,k)
        end do
       end do 
      end if
      if(itest.eq.2) then
       do i=1,ntrain
        do k=1,nft
         a4(i,k)=a(i,k)
        end do
       end do 
       nl=nrecord-ntrain
       do i=1,nl
        nlist(i)=i+ntrain
       end do
      end if
      nrecord2=nl
      nrecord1=nrecord-nrecord2
c=======================================================================
      if(nrecord2.gt.0) then
       knn=nrecord1/100
       knn=max(knn,3)
       knn=min(knn,17)
       do i=1,nl
        print *,i
        i1=nlist(i)
        d6=0.0d+32
        d7=0.0d+00
        do j=1,nrecord
         do i2=1,nl
          IF(j.eq.nlist(i2)) GO TO 51 
         end do
         d3(j)=0.0d+00
         do k=1,nft-1
          d3(j)=d3(j)+(a(i1,k)-a(j,k))**2
         end do
         d6=dmin1(d6,d3(j))
         d7=dmax1(d7,d3(j))         
  51    end do
        eps=0.0d+00
  61    d4=d6+eps*(d7-d6)
        j1=0
        do j=1,nrecord
         do i2=1,nl
          IF(j.eq.nlist(i2)) GO TO 71
         end do
         if(d3(j).le.d4) then
          j1=j1+1
          listtest(i,j1)=j
          if(j1.ge.knn) go to 81
         end if
  71    end do           
        eps=eps+1.0d-02
        go to 61 
  81   end do
      end if
c=======================================================================
      call scaling
c=======================================================================
      ncount=0
      tlin=0.0d+00
      do nc=1,nclust
       PRINT 3,nc
 3     FORMAT('         The number of hyperplanes:',i5)
       if (nc.eq.1) then
        do j=1,nft
         z(j)=1.0d+00
        end do
        nel1=nrecord1
        do k=1,nrecord1
         nob1(k)=k
        end do
        call optimbfgs(z,z1)
        call clusterfunc(f,z1)
        ncount=ncount+1
        do j=1,nft
         x(j)=z1(j)
        end do
        GO TO 2
       end if

       call step1(x,xc,ngood0)
       fmin=1.0d+32
       do i=1,ngood0
        do j=1,nft
         z(j)=xc(j+(i-1)*nft)
        end do

        nel1=0
        do k=1,nrecord1
         f2=z(nft)
         do j=1,nft-1
          f2=f2+a4(k,j)*z(j)
         end do
         f3=(a4(k,nft)-f2)**2
         if (f3.lt.dminim(k)) then
            nel1=nel1+1
            nob1(nel1)=k
         end if
        end do
        if (nel1.gt.0) then
         call optimbfgs(z,z1)
         call auxfunc(z1,f)
         ncount=ncount+1
         do j=1,nft
          xc(j+(i-1)*nft)=z1(j)
         end do
        end if
        if (nel1.eq.0) f=1.0d+32
        fv2(i)=f
        fmin=dmin1(fmin,f)
       end do

       jpoints=0
       d10=gamma2*fmin
       do i=1,ngood0
        fv1=fv2(i)
        if (fv2(i).le.d10) then
         jpoints=jpoints+1
         do j=1,nft
          xc(j+(jpoints-1)*nft)=xc(j+(i-1)*nft)
         end do
        end if
       end do
       call cleaning(jpoints,xc)
       ngood0=jpoints
       fb=1.0d+32
       do i=1,ngood0
        do j=1,nft
         x(j+(nc-1)*nft)=xc(j+(i-1)*nft)
        end do
        call optim(x,x1,f)
        ncount=ncount+nc
        if (f.lt.fb) then
         fb=f
         do k=1,nc
          do j=1,nft
           xbest(j+(k-1)*nft)=x1(j+(k-1)*nft)
          end do
         end do
        end if
       end do

       do k=1,nc
        do j=1,nft
         x(j+(k-1)*nft)=xbest(j+(k-1)*nft)
        end do
       end do

  2    call distribution(x,f)
       IF((nc.eq.1).and.(nrecord1.le.1000)) toler=1.0d-08*f
       IF((nc.eq.1).and.(nrecord1.gt.1000)) toler=5.0d-06*f
       call rescaling(x,x4)        

       call checkfinal(x4,ff,rmse1,rmse2,rmse3,rmse4,rmse5,rmae1
     1 ,rmae2,rmae3,rmae4,rmae5,ce1,ce2,ce3,ce4,ce5,pirson1,pirson2
     2 ,pirson3,pirson4,pirson5)
       
       ftr(nc,nfold)=ff
       
       rmsetr(nc,nfold)=rmse1
       rmsetest2(nc,nfold)=rmse2
       rmsetest3(nc,nfold)=rmse3
       rmsetest4(nc,nfold)=rmse4
       rmsetest5(nc,nfold)=rmse5
       
       rmaetr(nc,nfold)=rmae1
       rmaetest2(nc,nfold)=rmae2
       rmaetest3(nc,nfold)=rmae3
       rmaetest4(nc,nfold)=rmae4
       rmaetest5(nc,nfold)=rmae5
       
       cetr(nc,nfold)=ce1
       cetest2(nc,nfold)=ce2
       cetest3(nc,nfold)=ce3
       cetest4(nc,nfold)=ce4
       cetest5(nc,nfold)=ce5
       
       pirs1(nc,nfold)=pirson1
       pirs2(nc,nfold)=pirson2
       pirs3(nc,nfold)=pirson3
       pirs4(nc,nfold)=pirson4
       pirs5(nc,nfold)=pirson5
       
       write(40,*)
       write(40,573) nc
 573   format('No of hyperplanes:',i4)
       write(40,*)
       write(41,*)
       write(41,*)
       do k=1,nc
        WRITE(41,722) (x4(j+(k-1)*nft),j=1,nft)
       end do
722    FORMAT(20f12.5)
       write(40,543) ff
 543   format('Fit function value:',f24.8)
       WRITE(40,*)
       write(40,541) ncount
 541   FORMAT('The number of linear regression problems solved:',i10)
       tlin1=tlin/dble(nrecord1)
       WRITE(40,*)
       write(40,*)
       write(40,511) tlin
 511   FORMAT('The number of linear function evaluations:',f20.0)
       WRITE(40,*)
       write(40,519) tlin1
 519   FORMAT('Average number of linear function evaluations:',f20.0)
       WRITE(40,*)
       call cpu_time(time4)
       timef=time4-time1
       write(42,610) nc,ff,ncount,tlin,tlin1,timef
 610   FORMAT(I5,f26.6,I12,f28.0,f24.0,f13.3)
       WRITE(40,*)
       write(40,141) timef
       WRITE(40,*)
       write(40,*)
       WRITE(40,*)
       if(timef.gt.tlimit) go to 210
      end do
141   FORMAT('CPU time:',f12.3)
      GO TO 211
 210  CONTINUE
      write(43,803) 
 803  FORMAT('Fit function values on training set:')
      WRITE(43,*)
      do i=1,nclust
       write(43,332) i
       WRITE(43,331) (ftr(i,j),j=1,nfoldmax)
      end do
 331  FORMAT(5f28.6)
 332  FORMAT('#Clusters:',i6)
      WRITE(43,*) 
      write(43,804) 
 804  FORMAT('VALUES OF PERFORMANCE MEASURES ON TRAINING SET:')       
      WRITE(43,*)      
      do i=1,nclust
       WRITE(43,*)
       write(43,352) i
       write(43,*)
       WRITE(43,351) (rmsetr(i,j),j=1,nfoldmax)
       WRITE(43,354) (rmaetr(i,j),j=1,nfoldmax)
       WRITE(43,353) (cetr(i,j),j=1,nfoldmax)
       WRITE(43,355) (pirs1(i,j),j=1,nfoldmax)
      end do
 351  FORMAT('RMSE-tr:  ',5f28.6)
 354  FORMAT('MAE-tr:   ',5f28.6)
 353  FORMAT('CD-tr:    ',5f28.6)
 355  FORMAT('PIRSON-tr:',5f28.6)
 352  FORMAT('#Clusters:',i6) 
       
      if(itest.ge.2) then
       WRITE(43,*)
       WRITE(43,809)
 809   FORMAT('__________________________________________________')      
       WRITE(43,*)       
       write(43,805) 
 805   FORMAT('VALUES OF PERFORMANCE MEASURES ON TEST SET:')         
       WRITE(43,*)      
       do i=1,nclust
        WRITE(43,*)
        write(43,362) i
        write(43,*)
        write(43,806) 
 806    FORMAT('The first method, largest cluster:')         
        write(43,*)  
        WRITE(43,361) (rmsetest2(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest2(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest2(i,j),j=1,nfoldmax)
        WRITE(43,365) (pirs2(i,j),j=1,nfoldmax)
        write(43,*)        
        write(43,807) 
 807    FORMAT('The second method, using weights:')         
        write(43,*)  
        WRITE(43,361) (rmsetest3(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest3(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest3(i,j),j=1,nfoldmax)
        WRITE(43,365) (pirs3(i,j),j=1,nfoldmax)
        write(43,*)
        write(43,808) 
 808    FORMAT('The third method, k-NN:')         
        write(43,*)  
        WRITE(43,361) (rmsetest4(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest4(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest4(i,j),j=1,nfoldmax)
        WRITE(43,365) (pirs4(i,j),j=1,nfoldmax)
        write(43,*)
        write(43,810) 
 810    FORMAT('The fourth method, distances:')         
        write(43,*)  
        WRITE(43,361) (rmsetest5(i,j),j=1,nfoldmax)
        WRITE(43,364) (rmaetest5(i,j),j=1,nfoldmax)
        WRITE(43,363) (cetest5(i,j),j=1,nfoldmax)
        WRITE(43,365) (pirs5(i,j),j=1,nfoldmax)
       end do
 361   FORMAT('RMSE-test:  ',5f28.6)
 364   FORMAT('MAE-test:   ',5f28.6)
 363   FORMAT('CD-test:    ',5f28.6)
 365   FORMAT('PIRSON-test:',5f28.6)
 362   FORMAT('#Clusters:',i6)      
      end if

      close(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      close(78)
      stop
      end

c=======================================================================
c     distribution over clusters
c=======================================================================
      subroutine distribution(x,f)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),dminim(maxrec),a4(maxrec,maxnft)
      integer nob(maxclust,maxrec),nel(maxclust)
      common /cnft/nft,/cnk/nel,nob,/cdminim/dminim,/cnc/nc
     1 ,/crecord1/nrecord1,/c24/a4,/ctlin/tlin

      f=0.0d+00
      do i=1,nc
       nel(i)=0
      END do

      do k=1,nrecord1
       d1=1.0d+12
       do i=1,nc
        d2=x(i*nft)
        do j=1,nft-1
         d2=d2+a4(k,j)*x(j+(i-1)*nft)
        end do
        tlin=tlin+1.0d+00
        d3=(d2-a4(k,nft))**2
        if (d3.lt.d1) then
         d1=d3
         i1=i
        end if
       end do
       f=f+d1
       nel(i1)=nel(i1)+1
       nob(i1,nel(i1))=k
       dminim(k)=d1
      end do
c=======================================================================
      return
      end

c=======================================================================
c  final distribution over clusters
c=======================================================================
      subroutine checkfinal(x,ff,rmse1,rmse2,rmse3,rmse4,rmse5,rmae1
     1 ,rmae2,rmae3,rmae4,rmae5,ce1,ce2,ce3,ce4,ce5,pirson1,pirson2
     2 ,pirson3,pirson4,pirson5)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),w(maxclust)
     1 ,ebar(maxrec),ebar2(maxrec),xw(maxnft),xc(maxdim),dist(maxclust)
      integer nlist(maxrec),nob(maxclust,maxrec),nel(maxclust)
     1 ,listtest(100000,50)
      common /c22/a,/cnft/nft,/cnlist/nl,nlist,/cnrecord2/nrecord2
     1 ,/crecord1/nrecord1,/crecord/nrecord,/cnc/nc,/cnk/nel,nob
     2 ,/cltest/knn,listtest
c=======================================================================
      r1=0.0d+00
      do i=1,nrecord
       do i1=1,nl
        IF(i.eq.nlist(i1)) go to 1 
       end do
       r1=r1+a(i,nft)
  1   end do
  
      r2=0.0d+00  
      do i=1,nl
       i1=nlist(i)
       r2=r2+a(i1,nft)
      end do
 
      r1=r1/dble(nrecord1)
      if(nrecord2.gt.0) r2=r2/dble(nrecord2)
      
      fe1=0.0d+00
      do i=1,nrecord
       do i1=1,nl
        IF(i.eq.nlist(i1)) go to 2
       end do
       fe1=fe1+(a(i,nft)-r1)**2
  2   end do
  
      fe2=0.0d+00
      do i=1,nl
       i1=nlist(i) 
       fe2=fe2+(a(i1,nft)-r2)**2
      end do
c=======================================================================
c calculation of performance measures for training set
c=======================================================================
      ebar0=0.0d+00
      ff=0.0d+00    
      do i=1,nrecord
       do i1=1,nl
        IF(i.eq.nlist(i1)) GO TO 3
       end do
       f4=1.0d+26
       do k=1,nc
        f5=x(k*nft)
        do j=1,nft-1
         f5=f5+a(i,j)*x(j+(k-1)*nft)
        end do
        f3=(f5-a(i,nft))**2
        if(f3.lt.f4) then 
         f4=f3
         ebar(i)=f5
        end if
       end do
       ff=ff+f4
       ebar0=ebar0+ebar(i)
  3   end do
      ebar0=ebar0/dble(nrecord1)

      f=0.0d+00
      rmae1=0.0d+00
      r3=0.0d+00
      r4=0.0d+00
      do i=1,nrecord
       do i1=1,nl
        IF(i.eq.nlist(i1)) GO TO 4 
       end do
       f=f+(ebar(i)-a(i,nft))**2
       rmae1=rmae1+dabs(ebar(i)-a(i,nft))
       r3=r3+(a(i,nft)-r1)*(ebar(i)-ebar0)
       r4=r4+(ebar(i)-ebar0)**2
  4   end do
      rmse1=sqrt(f/DBLE(nrecord1))
      rmae1=rmae1/DBLE(nrecord1)
      pirson1=r3/sqrt(fe1*r4)
      ce1=1.0d+00-f/fe1
c=======================================================================
c cal. of per. measures for test set - the use of the largest cluster 
c=======================================================================
      if(nrecord2.gt.0) then
       k1=0
       do i=1,nc
        if(k1.lt.nel(i)) then
         k1=nel(i)
         kmax=i
        end if
       end do

       ebar20=0.0d+00    
       do i=1,nl
        i1=nlist(i)
        d5=x(kmax*nft)
        do j=1,nft-1
         d5=d5+a(i1,j)*x(j+(kmax-1)*nft)
        end do
        ebar2(i)=d5
        ebar20=ebar20+ebar2(i)
       end do      
       ebar20=ebar20/dble(nl)

       ftest=0.0d+00
       rmae2=0.0d+00
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)
        ftest=ftest+(ebar2(i)-a(i1,nft))**2
        rmae2=rmae2+dabs(ebar2(i)-a(i1,nft))
        r3=r3+(a(i1,nft)-r2)*(ebar2(i)-ebar20)
        r4=r4+(ebar2(i)-ebar20)**2
       end do
       rmse2=sqrt(ftest/DBLE(nl))
       rmae2=rmae2/DBLE(nl)
       pirson2=r3/sqrt(fe2*r4)
       ce2=1.0d+00-ftest/fe2
      end if 
c=======================================================================
c cal. of per. measures for test set - the use of weights 
c=======================================================================
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

       ebar20=0.0d+00    
       do i=1,nl
        i1=nlist(i)
        d5=xw(nft)
        do j=1,nft-1
         d5=d5+a(i1,j)*xw(j)
        end do
        ebar2(i)=d5
        ebar20=ebar20+ebar2(i)
       end do      
       ebar20=ebar20/dble(nl)

       ftest=0.0d+00
       rmae3=0.0d+00
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)
        ftest=ftest+(ebar2(i)-a(i1,nft))**2
        rmae3=rmae3+dabs(ebar2(i)-a(i1,nft))
        r3=r3+(a(i1,nft)-r2)*(ebar2(i)-ebar20)
        r4=r4+(ebar2(i)-ebar20)**2
       end do
       rmse3=sqrt(ftest/DBLE(nl))
       rmae3=rmae3/DBLE(nl)
       pirson3=r3/sqrt(fe2*r4)
       ce3=1.0d+00-ftest/fe2
      end if 
c======================================================================
c cal. of per. measures for test set - the use of k nearest neighbors 
c======================================================================
      if(nrecord2.gt.0) then

       do i=1,nl
        i1=nlist(i)
        do k=1,nc
         w(k)=0.0d+00
        end do
        
        do j=1,knn
         j1=listtest(i,j)
         do k=1,nc
          do k1=1,nel(k)
           k2=nob(k,k1)
           if(j1.eq.k2) then
            w(k)=w(k)+1.0d+00/dble(knn)
            go to 9
           end if
          end do
         end do
   9    end do
       
        do j=1,nft
         xw(j)=0.0d+00
         do k=1,nc
          xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
         end do
        end do

        d5=xw(nft)
        do j=1,nft-1
         d5=d5+a(i1,j)*xw(j)
        end do
        ebar2(i)=d5
       end do

       ftest=0.0d+00
       rmae4=0.0d+00
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)
        ftest=ftest+(ebar2(i)-a(i1,nft))**2
        rmae4=rmae4+dabs(ebar2(i)-a(i1,nft))
        r3=r3+(a(i1,nft)-r2)*(ebar2(i)-ebar20)
        r4=r4+(ebar2(i)-ebar20)**2
       end do
       rmse4=sqrt(ftest/DBLE(nl))
       rmae4=rmae4/DBLE(nl)
       pirson4=r3/sqrt(fe2*r4)
       ce4=1.0d+00-ftest/fe2
      end if
c======================================================================
c cal. of per. measures for test set - the use of distances 
c======================================================================
      if(nrecord2.gt.0) then
       do k=1,nc
        do j=1,nft
         xc(j+(k-1)*nft)=0.0d+00
        end do
       end do
       do k=1,nc
        k1=nel(k)
        do j=1,k1
         k2=nob(k,j)
         do j1=1,nft-1
          xc(j1+(k-1)*nft)=xc(j1+(k-1)*nft)+a(k2,j1)
         end do
        end do
       end do
       do k=1,nc
        if(nel(k).gt.0) then
         do j=1,nft
          xc(j+(k-1)*nft)=xc(j+(k-1)*nft)/dble(nel(k))
         end do
        end if
       end do      

       do i=1,nl
        i1=nlist(i)
        do k=1,nc
         dist(k)=0.0d+00
         do j=1,nft-1
          dist(k)=dist(k)+(a(i1,j)-xc(j+(k-1)*nft))**2
         end do
         dist(k)=sqrt(dist(k))
        end do
        
        d0=0.0d+00
        do k=1,nc
         d0=d0+1.0d+00/dist(k)**2
        end do

        do k=1,nc
         w(k)=1.0d+00/(dist(k)**2*d0)
        end do        

        do j=1,nft
         xw(j)=0.0d+00
         do k=1,nc
          xw(j)=xw(j)+w(k)*x(j+(k-1)*nft)
         end do
        end do
       
        d5=xw(nft)
        do j=1,nft-1
         d5=d5+a(i1,j)*xw(j)
        end do
        ebar2(i)=d5       
       end do 

       ftest=0.0d+00
       rmae5=0.0d+00
       r3=0.0d+00
       r4=0.0d+00
       do i=1,nl
        i1=nlist(i)
        ftest=ftest+(ebar2(i)-a(i1,nft))**2
        rmae5=rmae5+dabs(ebar2(i)-a(i1,nft))
        r3=r3+(a(i1,nft)-r2)*(ebar2(i)-ebar20)
        r4=r4+(ebar2(i)-ebar20)**2
       end do
       rmse5=sqrt(ftest/DBLE(nl))
       rmae5=rmae5/DBLE(nl)
       pirson5=r3/sqrt(fe2*r4)
       ce5=1.0d+00-ftest/fe2
      end if
c=======================================================================
      return
      end

c=======================================================================
c Step1 initialize clusters
c=======================================================================
      subroutine step1(x,z,jpoints)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100
     1 ,maxinit=1000000)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a4(maxrec,maxnft),z(maxinit)
     1 ,dminim(maxrec),d1(maxclust,maxrec),d5(maxrec)
      INTEGER nel(maxclust),nob(maxclust,maxrec),l1(maxrec),l2(maxrec)
     1 ,l3(maxrec) 
      common /c24/a4,/cnft/nft,/crecord1/nrecord1,/cnk/nel,nob
     1 ,/cgamma/gamma1,/cdminim/dminim,/cnc/nc,/ctoler/toler

      eps=5.0d+00*toler

      do i=1,nc-1
       do k=1,nel(i)
        k1=nob(i,k)
        l1(k1)=i
       end do
      end do

      do i=1,nc-1
       do k=1,nrecord1
        d0=0.0d+00
        do j=1,nft-1
         d0=d0+x(j+(i-1)*nft)*a4(k,j)
        end do
        d1(i,k)=d0
       end do   
      end do
      
      l0=0
      do i=1,nc-1
       lnc=0
       do k=1,nel(i)
        k1=nob(i,k)
        do j=1,lnc
         k2=l3(j)
         d3=d1(i,k1)+a4(k2,nft)-d1(i,k2)
         d3=abs(d3-a4(k1,nft))
         if(d3.le.eps) go to 1
        end do
        lnc=lnc+1
        l3(lnc)=k1
  1    end do
       do k=1,lnc
        l2(k+l0)=l3(k)
       end do
       l0=l0+lnc
      end do

      d6=0.0d+00
      do i=1,l0
       i1=l2(i)
       i2=l1(i1)
       d5(i)=0.0d+00
       do k=1,nrecord1
        d2=d1(i2,k)+a4(i1,nft)-d1(i2,i1)
        d3=(d2-a4(k,nft))**2
        d4=dmin1(0.0d+00,d3-dminim(k))
        d5(i)=d5(i)+d4
       end do
       d6=dmin1(d6,d5(i))
      end do

      jpoints=0
      d3=gamma1*d6
      do k=1,l0
       i1=l2(k)
       i2=l1(i1)
       if (d5(k).le.d3) then
        jpoints=jpoints+1
        do j=1,nft-1
         z(j+(jpoints-1)*nft)=x(j+(i2-1)*nft)
        end do
        z(jpoints*nft)=a4(i1,nft)-d1(i2,i1)
       end if
      end do
      call cleaning(jpoints,z)
      return
      end

c=======================================================================
      subroutine cleaning(jpoints,z)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100
     1 ,maxinit=1000000)
      implicit double precision (a-h,o-z)
      double precision z(maxinit),z1(maxnft),z2(maxinit)
      common /cnft/nft,/ctoler/toler

      jpoints0=0
      do i=1,jpoints
       do j=1,nft
        z1(j)=z(j+(i-1)*nft)
       end do
       if (i.eq.1) then
        jpoints0=1
        do j=1,nft
         z2(j)=z1(j)
        end do
       end if

       if (i.gt.1) then
        do k=1,jpoints0
         d2=0.0d+00
         do j=1,nft
          d2=d2+ABS(z1(j)-z2(j+(k-1)*nft))
         end do
         IF(d2.le.toler) GO TO 1
        end do
        jpoints0=jpoints0+1
        do j=1,nft
         z2(j+(jpoints0-1)*nft)=z1(j)
        end do
       end if
   1  end do
      do k=1,jpoints0
       do j=1,nft
        z(j+(k-1)*nft)=z2(j+(k-1)*nft)
       end do
      end do
      jpoints=jpoints0
      end subroutine

c======================================================================
      subroutine scaling
      PARAMETER(maxnft=100, maxrec=500000)
      implicit double precision(a-h,o-z)
      double precision a4(maxrec,maxnft),clin(maxnft),dlin(maxnft)
      common /c24/a4,/crecord1/nrecord1,/cnft/nft,/ctransform/clin,dlin
c======================================================================
      do i=1,nft
       dm=0.0d+00
       do j=1,nrecord1
        dm=dm+a4(j,i)
       END do
       dm=dm/DBLE(nrecord1)
       var=0.0d+00
       do j=1,nrecord1
         var=var+(a4(j,i)-dm)**2
       end do
       var=var/DBLE(nrecord1-1)
       var=dsqrt(var)
       if(var.ge.1.0d-04) then
            clin(i)=1.0d+00/var
            dlin(i)=-dm/var
            do j=1,nrecord1
              a4(j,i)=clin(i)*a4(j,i)+dlin(i)
            end do
       end if
       if (var.lt.1.0d-04) then
         if (dabs(dm).le.1.0d+00) then
          clin(i)=1.0d+00
          dlin(i)=0.0d+00
         end if
         if (dabs(dm).gt.1.0d+00) then
          clin(i)=1.0d+00/dm
          dlin(i)=0.0d+00
         end if
         do j=1,nrecord1
          a4(j,i)=clin(i)*a4(j,i)+dlin(i)
         end do
       end if
      end do
c=======================================================================
      return
      end

c=======================================================================
c  rescaling
c=======================================================================
      subroutine rescaling(x,x4)
      PARAMETER(maxdim=5000, maxnft=100, maxrec=500000)
      implicit double precision(a-h,o-z)
      double precision clin(maxnft),dlin(maxnft),x(maxdim),x4(maxdim)
      common /cnft/nft,/ctransform/clin,dlin,/cnc/nc
c=======================================================================
      rabs=dabs(clin(nft))
      if (rabs.gt.1.0d-08) then
       do i=1,nc
        d0=0.0d+00
        do j=1,nft-1
         x4(j+(i-1)*nft)=x(j+(i-1)*nft)*clin(j)/rabs
         d0=d0+x(j+(i-1)*nft)*dlin(j)
        end do
        x4(i*nft)=(x(i*nft)+d0-dlin(nft))/rabs
       end do
      end if
      return
      end

c=======================================================================
      subroutine auxfunc(x,fval)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
      double precision x(maxdim),a4(maxrec,maxnft),dminim(maxrec)
      common /c24/a4,/crecord1/nrecord1,/cnft/nft,/cdminim/dminim
     1 ,/ctlin/tlin
      fval=0.0d+00
      do i=1,nrecord1
       f1=dminim(i)
       f2=x(nft)
       do j=1,nft-1
        f2=f2+a4(i,j)*x(j)
       end do
       tlin=tlin+1.0d+00
       f3=(a4(i,nft)-f2)**2
       fval=fval+dmin1(f1,f3)
      end do
      return
      end

c=======================================================================
      subroutine clusterfunc(f,x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
      double precision x(maxdim),a4(maxrec,maxnft)
      common /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/ctlin/tlin
      f=0.0d+00
      do i=1,nrecord1
       f12=1.0d+32
       do k=1,nc
        d1=x(k*nft)
        do j=1,nft-1
         d1=d1+a4(i,j)*x(j+(k-1)*nft)
        end do
        tlin=tlin+1.0d+00
        d2=(d1-a4(i,nft))**2
        f12=dmin1(f12,d2)
       end do
       f=f+f12
      end do
      return
      end

c=======================================================================
      subroutine clusterfunc1(f,x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
      double precision x(maxdim),a4(maxrec,maxnft)
      common /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/ctlin/tlin
      f=0.0d+00
      do i=1,nrecord1
       do k=1,nc
        d1=x(k*nft)
        do j=1,nft-1
         d1=d1+a4(i,j)*x(j+(k-1)*nft)
        end do
        tlin=tlin+1.0d+00
        d2=(d1-a4(i,nft))**2
        f=f+d2
       end do
      end do
      return
      end
      
      subroutine clusterfunc2(f,x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
      double precision x(maxdim),a4(maxrec,maxnft),d2(maxclust)
      common /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/ctlin/tlin
      f=0.0d+00
      do i=1,nrecord1
        d5=0.0d+00
        do k=1,nc
         d1=x(nft*k)
         do j=1,nft-1
          d1=d1+x(j+(k-1)*nft)*a4(i,j)
         end do
         tlin=tlin+1.0d+00
         d2(k)=(d1-a4(i,nft))**2
         d5=d5+d2(k)
        end do      
        d4=0.0d+00
        do k=1,nc
         d6=d5-d2(k)
         if (d4.lt.d6) d4=d6
        end do      
       f=f+d4
      end do
      return
      end      
      
c=======================================================================
      subroutine fgrad1(x1,grad1)
       implicit double precision (a-h,o-z)
       PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
       double precision a4(maxrec,maxnft), grad1(maxdim), x1(maxdim)
      common /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/csize/m
     1 ,/ctlin/tlin

      do i=1,m
       grad1(i)=0.0d+00
      end do

      do i=1,nrecord1
       do k=1,nc
        d1=x1(k*nft)
        do j=1,nft-1
         d1=d1+a4(i,j)*x1(j+(k-1)*nft)
        end do
        tlin=tlin+1.0d+00
        d1=d1-a4(i,nft)
        do j=1,nft-1
         grad1(j+(k-1)*nft)=grad1(j+(k-1)*nft)+2.0d+00*d1*a4(i,j)
        end do
        grad1(k*nft)=grad1(k*nft)+2.0d+00*d1
       end do
      END do
      end subroutine

      subroutine fgrad2(x,grad2)
       implicit double precision (a-h,o-z)
       PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
       double precision x(maxdim), a4(maxrec,maxnft), d2(maxclust)
     1  ,dminim(maxrec),grad2(maxdim),d3(maxclust)
      common /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/csize/m
     1 ,/cdminim/dminim,/ctlin/tlin
     
      do i=1,m
       grad2(i)=0.0d+00
      end do

      do i=1,nrecord1
       do k=1,nc
        d3(k)=x(k*nft)
        do j=1,nft-1
         d3(k)=d3(k)+a4(i,j)*x(j+(k-1)*nft)
        end do
        tlin=tlin+1.0d+00
        d2(k)=(d3(k)-a4(i,nft))**2
       end do
       d4=0.0d+00
       do j=1,nc
        d5=0.0d+00
        do k=1,nc
         if (k.ne.j) then
          d5=d5+d2(k)
         end if
        end do
        if (d4.lt.d5) then
         d4=d5
         jindex=j
        end if
       end do

       do j=1,nc
        if (j.ne.jindex) then
         do k=1,nft-1
          grad2(k+(j-1)*nft)=grad2(k+(j-1)*nft)
     1       +2.0d+00*(d3(j)-a4(i,nft))*a4(i,k)
         end do
         grad2(j*nft)=grad2(j*nft)+2.0d+00*(d3(j)-a4(i,nft))
        end if
       end do
      end do
      end subroutine

      SUBROUTINE optim(xbar,x,fvalue)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000)
      double precision x0(maxvar),x(maxvar),grad1(maxvar),grad2(maxvar)
     1 ,gradbar(maxvar),xbar(maxvar) 
      COMMON /cnft/nft,/cnc/nc,/cnopt/nopt,/cgbar/gradbar,/cfbar/fbar
     1 ,/csize/m
c=======================================================================
c  Input data:
c  n        - number of variables
c=======================================================================
      n = nc*nft
      eps=1.0d-02
      eps2=1.0d-05
      m=n
      nopt=1
      iter=0
  1   continue
      iter=iter+1
      call fgrad1(xbar,grad1)
      call fgrad2(xbar,grad2)
      d1=0.0d+00
      do k=1,m
       d1=d1+(grad1(k)-grad2(k))**2
      end do
      d1=sqrt(d1)       
      if(d1.le.eps) go to 3
      call clusterfunc(fvalue,xbar)
      if(iter.gt.1) then
       dif=abs(fold-fvalue)/(fold+1.0d+00)
       if(dif.le.eps2) go to 3
      end if      
      fold=fvalue
      do j=1,nc
       do k=1,nft
        x0(k)=xbar(k+(j-1)*nft)
        gradbar(k)=grad2(k+(j-1)*nft) 
       end do
       d2=0.0d+00
       do k=1,nft
        d2=d2+gradbar(k)*x0(k)
       end do
       fbar=-d2
       call bfgs(nft,x0,x)
       do k=1,nft
        xbar(k+(j-1)*nft)=x(k)
       end do
      end do
      go to 1
 3    continue
      do j=1,m
       x(j)=xbar(j)
      end do
      call clusterfunc(fvalue,x)
c=======================================================
      RETURN
      END

c=======================================================================
      SUBROUTINE optimbfgs(x0,x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000)
      double precision x0(maxvar),x(maxvar)
      COMMON /cnft/nft,/cnopt/nopt,/csize/m
c=======================================================================
c  Input data:
c  n        - number of variables
c=======================================================================
      n=nft
      m=n
      nopt=2
c=======================================================================
c  Calling BFGS
c=======================================================================
      call bfgs(n,x0,x)
c=======================================================================
      RETURN
      END

c=======================================================================
      DOUBLE PRECISION FUNCTION func(x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000, maxrec=500000, maxnft=100)
      double precision x(maxvar),a4(maxrec,maxnft),grad2(maxvar)
      INTEGER nob1(maxrec)
      COMMON /cnel1/nel1,nob1,/c24/a4,/cnft/nft,/cnumreg/tnfreg,tngreg
     1 ,/cnc/nc,/cnopt/nopt,/cfbar/fbar,/cgbar/grad2,/ctlin/tlin
c=======================================================================
      if(nopt.eq.1) then
       call clusterfunc3(f1,x)
       f2=0.0d+00
       do j=1,nft
        f2=f2+grad2(j)*x(j)
       end do
       f2=f2+fbar
       func=f1-f2
      end if
      
      if(nopt.eq.2) then
       d2=0.0d+00
       do k=1,nel1
        d1=x(nft)
        do j=1,nft-1
         d1=d1+a4(nob1(k),j)*x(j)
        end do
        tlin=tlin+1.0d+00
        d2=d2+(d1-a4(nob1(k),nft))**2
       end do
       func=d2
       tnfreg=tnfreg+dble(nel1)
      end if
c=======================================================================
      return
      end

c=======================================================================
      subroutine clusterfunc3(f,x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
      double precision x(maxdim),a4(maxrec,maxnft)
      common /c24/a4,/crecord1/nrecord1,/cnc/nc,/cnft/nft,/ctlin/tlin
      f=0.0d+00
      do i=1,nrecord1
       d1=x(nft)
       do j=1,nft-1
        d1=d1+a4(i,j)*x(j)
       end do
       tlin=tlin+1.0d+00
       d2=(d1-a4(i,nft))**2
       f=f+d2
      end do
      return
      end

c======================================================================
      subroutine fgrad3(x1,grad1)
       implicit double precision (a-h,o-z)
       PARAMETER(maxdim=5000, maxrec=500000, maxclust=50, maxnft=100)
       double precision a4(maxrec,maxnft), grad1(maxdim), x1(maxdim)
      common /c24/a4,/crecord1/nrecord1,/cnft/nft,/ctlin/tlin

      do i=1,nft
       grad1(i)=0.0d+00
      end do

      do i=1,nrecord1
       d1=x1(nft)
       do j=1,nft-1
        d1=d1+a4(i,j)*x1(j)
       end do
       tlin=tlin+1.0d+00
       d1=d1-a4(i,nft)
       do j=1,nft-1
        grad1(j)=grad1(j)+2.0d+00*d1*a4(i,j)
       end do
       grad1(nft)=grad1(nft)+2.0d+00*d1
      END do
      end subroutine

c=======================================================================
      subroutine dfunc(x,grad)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000, maxrec=500000, maxclust=50, maxnft=100)
      double precision x(maxvar),grad(maxvar),a4(maxrec,maxnft)
     1 ,grad1(maxvar),grad2(maxvar) 
      INTEGER nob1(maxrec)
      COMMON /cnel1/nel1,nob1,/c24/a4,/cnft/nft,/cnumreg/tnfreg,tngreg
     1 ,/cnc/nc,/cnopt/nopt,/cgbar/grad2,/ctlin/tlin
c=======================================================================
      if(nopt.eq.1) then
       call fgrad3(x,grad1)
       do j=1,nft
        grad(j)=grad1(j)-grad2(j)
       end do
      end if

      if(nopt.eq.2) then
       do j=1,nft
        grad(j)=0.0d+00
       end do
       do k=1,nel1
        d1=x(nft)
        do j=1,nft-1
         d1=d1+a4(nob1(k),j)*x(j)
        end do
        tlin=tlin+1.0d+00
        d1=d1-a4(nob1(k),nft)
        do j=1,nft-1
         grad(j)=grad(j)+2.0d+00*d1*a4(nob1(k),j)
        end do
        grad(nft)=grad(nft)+2.0d+00*d1
       end do
       tngreg=tngreg+dble(nel1)
      end if
c=======================================================================
      return
      end

c=======================================================================
      subroutine bfgs(nvar,x0,p)
       implicit double precision (a-h,o-z)
       PARAMETER(maxvar=5000)
       double precision x0(maxvar),p(maxvar)
       EXTERNAL func,dfunc
       ndim=nvar
       GTOL=1.0d-03
       do i=1,ndim
        p(i)=x0(i)
       end do 
       call dfpmin(p,NDIM,GTOL,iter,fret,func,dfunc)
      return
      END
      
      SUBROUTINE dfpmin(p,n,gtol,iter,fret,func,dfunc)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=5000,ITMAX=1000,STPMX=1.0d+02,
     * EPS=5.0d-07,TOLX=4.0d+00*EPS)
      DOUBLE PRECISION p(n)
      DOUBLE PRECISION dg(NMAX),g(NMAX),hdg(NMAX),hessin(NMAX,NMAX)
     * ,pnew(NMAX),xi(NMAX)
      EXTERNAL dfunc,func
CU    USES dfunc,func,lnsrch
      LOGICAL check
      fp=func(p)
      call dfunc(p,g)
      sum1=0.0d+00
      do i=1,n
        do j=1,n
          hessin(i,j)=0.0d+00
        end do
        hessin(i,i)=1.0d+00
        xi(i)=-g(i)
        sum1=sum1+p(i)**2
      end do
      stpmax=STPMX*dmax1(sqrt(sum1),dble(n))
      do its=1,ITMAX
        iter=its
        call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,func)
        fp=fret
        do i=1,n
          xi(i)=pnew(i)-p(i)
          p(i)=pnew(i)
        end do
        test=0.0d+00
        do i=1,n
          temp=dabs(xi(i))/dmax1(dabs(p(i)),1.0d+00)
          if(temp.gt.test)test=temp
        end do
        if(test.lt.TOLX)return
        do i=1,n
          dg(i)=g(i)
        end do
        call dfunc(p,g)
        test=0.0d+00
        den=dmax1(fret,1.0d+00)
        do i=1,n
          temp=dabs(g(i))*dmax1(dabs(p(i)),1.0d+00)/den
          if(temp.gt.test)test=temp
        end do
        if(test.lt.gtol)return
        do i=1,n
          dg(i)=g(i)-dg(i)
        end do
        do i=1,n
          hdg(i)=0.0d+00
          do j=1,n
            hdg(i)=hdg(i)+hessin(i,j)*dg(j)
          end do
        end do
        fac=0.0d+00
        fae=0.0d+00
        sumdg=0.0d+00
        sumxi=0.0d+00
        do i=1,n
          fac=fac+dg(i)*xi(i)
          fae=fae+dg(i)*hdg(i)
          sumdg=sumdg+dg(i)**2
          sumxi=sumxi+xi(i)**2
        end do
        if(fac**2.gt.EPS*sumdg*sumxi)then
          fac=1.0d+00/fac
          fad=1.0d+00/fae
          do i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
          end do
          do i=1,n
            do j=1,n
              hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+
     *fae*dg(i)*dg(j)
            end do
          end do
        endif
        do i=1,n
          xi(i)=0.0d+00
          do j=1,n
            xi(i)=xi(i)-hessin(i,j)*g(j)
          end do
        end do
      end do
      return
      END

      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      implicit double precision (a-h,o-z)
      LOGICAL check
      DOUBLE PRECISION g(n),p(n),x(n),xold(n)
      PARAMETER (ALF=1.0d-06,TOLX=1.0d-06)
      EXTERNAL func
CU    USES func
      check=.false.
      sum1=0.0d+00
      do i=1,n
        sum1=sum1+p(i)*p(i)
      end do
      sum1=dsqrt(sum1)
      if(sum1.gt.stpmax)then
        do i=1,n
          p(i)=p(i)*stpmax/sum1
        end do
      endif
      slope=0.0d+00
      do i=1,n
        slope=slope+g(i)*p(i)
      end do
      test=0.0d+00
      do i=1,n
        temp=dabs(p(i))/dmax1(dabs(xold(i)),1.0d+00)
        if(temp.gt.test)test=temp
      end do
      alamin=TOLX/test
      alam=1.0d+00
1     continue
        do i=1,n
          x(i)=xold(i)+alam*p(i)
        end do
        f=func(x)
        if(alam.lt.alamin)then
          do i=1,n
            x(i)=xold(i)
          end do
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.0d+00*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.0d+00)then
              tmplam=-slope/(2.0d+00*b)
            else
              disc=b*b-3.0d+00*a*slope
              tmplam=(-b+sqrt(disc))/(3.0d+00*a)
            endif
            if(tmplam.gt.5.0d-01*alam) tmplam=5.0d-01*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=dmax1(tmplam,1.0d-01*alam)
      goto 1
      END
