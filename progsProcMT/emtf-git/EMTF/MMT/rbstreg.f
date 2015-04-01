      subroutine rbstreg(x,y,r,nd,np,nch,b,bsv,work,cov,xxinv,iter_max,
     &    icvg,nomit)

 
ccc    general complex robust multiple regression; uses qr
ccc    ==> nd,np,nch = length of data vectors, # of predictors, # of data channels to fit
ccc    ==> x(nd,np) = complex predictor variables
ccc    ==> y(nd,nch) = complex data to fit (currently: regressed one at a time)
ccc    ==> nomit = number of data points omitted by multiplying by zero
ccc    <== r(nd,nch) = complex residuals
ccc    <== b(np,nch) = complex regression parameters   
ccc    <== cov(*) = covariance estimates for noise in predicted channels

      integer nd,np,nch,ldx,lwork,iw1,iw2,info,icvg
      complex x(nd,np),y(nd,nch),r(nd,nch),b(np,nch),
     &   bsv(np,nch)
      character*1 side,trans,diag,uplo
      complex cov(*),cdotc,xxinv(np,np)
      logical ldone
      real work(*)

      parameter (r0 = 1.4, tol=.001)

ccc       cfac is "correction factor" for Huber style sd estimate
ccc          Huber penalty, complex data ....
      cfac = 1./(2*(1.-(1.+r0)*exp(-r0) ))

ccc   compute complex q-r decompostion : using LAPACK
      nt = nd*nch
      nt2 = 2*nt
      np2 = 2*np
      nd2 = 2*nd
      ldx = nd
      lwork = 100*np 
ccc       work(1:np) = scalar factors of elementary reflectors (save)
ccc       work(np+1:np+lwork)
      iw1 = 1
      iw2 = 2*np+1
      iw3 = iw2 + 2*lwork
      call cgeqrf(nd,np,x,ldx,work(iw1),work(iw2),lwork,info)
c      write(*,*) 'info, optimal block size',info,work(iw2)

      do iter = 0,iter_max
ccc      put "corrected" data into array r
         call datawt(iter,nd,nch,y,r,r,work(iw3),r0)

ccc      compute Q^Y (uses full unitary matrix Q, data vectors (Y)
ccc             are stored in R)
         side = 'L'
         trans = 'C'
         ldr = ldx
         call cunmqr(side,trans,nd,nch,np,x,ldx,work(iw1),r,ldr
     &    ,work(iw2),lwork,info)

ccc     save old estimates, move Q^Y into b
ccc     blank upper part of Q^Y (leaving rotated residuals)
         do i = 1,nch
            do j = 1,np
               bsv(j,i) = b(j,i)
               b(j,i) = r(j,i)
               r(j,i) = (0.,0.)
            enddo
         enddo

ccc      compute residuals for corrected data
         side = 'L'
         trans = 'N'
         call cunmqr(side,trans,nd,nch,np,x,ldx,work(iw1),r,ldr
     &       ,work(iw2),lwork,info)

ccc      calculate error covariance
         ij = 0
         ndf = max(1,nd-nomit-np)
         do i = 1,nch
            do j = 1,i
               ij = ij + 1
               if(iter.gt.0) then
                  cov(ij) = cfac*cdotc(nd,r(1,j),1,r(1,i),1)/ndf
               else
                  cov(ij) = cdotc(nd,r(1,j),1,r(1,i),1)/ndf
               endif
            enddo
         enddo

ccc      move Q^Y back into R; zero out bottom part for computation
ccc      of predicted data
         do i = 1,nch
            do j = 1,np
               r(j,i) = b(j,i)
            enddo
            do j = np+1,nd
               r(j,i) = (0.,0.)
             enddo
          enddo

ccc      now compute predicted data (uses full unitary matrix Q)
         side = 'L'
         trans = 'N'
         call cunmqr(side,trans,nd,nch,np,x,ldx,work(iw1),r,ldr
     &       ,work(iw2),lwork,info)

ccc      compute new estimates B
         trans = 'N'
         diag = 'N'
         uplo = 'U'
         call ctrtrs(uplo,trans,diag,np,nch,x,ldx,b,np,info)

ccc      check to see if iterative scheme has converged
         if(iter.gt.0) then
            chngmx = 0.0
            ldone = .true.
            do i = 1,nch
               chng = snrm2(np2,bsv(1,i),1)
               do j = 1,np
                  bsv(j,i) = b(j,i)-bsv(j,i)
               enddo
               chng = snrm2(np2,bsv(1,i),1)/chng
               chngmx = max(chng,chngmx)
               ldone = ldone .and. ( chng.le.tol )
            enddo
         else if(iter_max.eq.0) then
            ldone = .true.
         else
            ldone = .false.
         endif

         if(ldone) then
            icvg = iter
ccc         compute XXinv
            do i = 1,np
               do j = 1,np
                  if(i.eq.j) then
                     xxinv(i,j) = (1.,0.)
                  else
                     xxinv(i,j) = (0.,0.)
                  endif
               enddo
            enddo
            trans = 'C'
            diag = 'N'
            uplo = 'U'
            call ctrtrs(uplo,trans,diag,np,np,x,ldx,xxinv,np,info)
            trans = 'N'
            diag = 'N'
            uplo = 'U'
            call ctrtrs(uplo,trans,diag,np,np,x,ldx,xxinv,np,info)
            return
         else
            do i = 1,nch
               ii = (i*(i+1))/2
               iw = iw3+i-1
               work(iw) = real(cov(ii))
            enddo
         endif
      enddo ! iter

ccc   if you get here, the iterative minimization didn't converge
      icvg = - nint(chngmx/tol)
      return
      end
c__________________________________________________________________
      subroutine datawt(iter,nd,nch,y,r,c,var,r0)

ccc   given data (y) and predicted (r), compute modified data and return in c
ccc    (c can overwrite r)
      integer iter,nd,nch
      complex y(nd,nch),r(nd,nch),c(nd,nch)
      real r0,var(nch)

      if(iter.eq.0) then
ccc       just copy data into r
         nt = nd*nch
         call ccopy(nt,y,1,c,1)
      else if(iter.gt.0) then
ccc     iterative minimization of Huber penalty functional
         do i = 1,nch
            r0s = sqrt(var(i))*r0
            do j = 1,nd
               if(abs(y(j,i)-r(j,i)).gt.r0s) then
                  w = r0s/abs(y(j,i)-r(j,i))
                  c(j,i) = w*y(j,i) + (1.-w)*r(j,i)
               else
                  c(j,i) = y(j,i)
               endif
            enddo
         enddo
      endif
      return
      end
