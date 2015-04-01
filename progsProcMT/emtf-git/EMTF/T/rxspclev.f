      subroutine rxspclev(xx,nf,nchstk,s,nu,ru,tol,lrbst,w,
     & nchin,lref,iter)

      parameter (nchmx_loc = 100)
c       robust cross spectral density stacking routine;
c      this routine uses the first two channels as reference
c       channels; the remaining nchstk-2 channels are the predicted outputs
c       in this version the weighting parameter for weighting by power
c       is an input arqument (wtpram)
c       nchin is the total number of data channels in all files, which
c        can exceed nchstk in the case of remote reference data

c     new version 1-26-88; substantial changes made to enhance protection
c     against extreme leverage points.
c     parameters:
c         varcor : needed to adjust scale estimate to allow for downweighted data
c                 points; this is just the expectation of the weighted residuals
c                 when the errors are iid unit variance gaussian and the huber loss
c                 function is used with r sub 0 = 1.5
c         edfl1,alpha,c1,c2: parameters needed to control downweighting of leverage
c                 points;  the following approach to limiting the effect of
c                 extreme leverage points is used: (1) a robust estimate of the
c                 horizontal magnetic field sdm is computed by ...
c                     more latter
c      nchstk:   num. of working channals=(nch+2) for remote reference case
c      nchin :   total channels of data read in
c      lref  :   logical parameter for different mode;
c                lref = 'true' -  remote reference
c                lref = 'false' - no remote reference
c      parameter(varcor=.7784,edfl1=1000.0,alpha=.5,c1=1000.0,c2=1000.0,
c     &              c3=1000.0,itmax = 20)
      parameter(varcor=.7784,edfl1=10.0,alpha=.5,c1=2.0,c2=20.0,
     &              c3=10.0,itmax = 50)
 
c       external functions wt,rdcndwt are weighting functions for
c       pulling residuals in
 
      external wt,rdcndwt
 
      integer sindex
      logical ldone(nchmx_loc),lastit,lrbst,lref
      real ru(*),scale(nchmx_loc),rus(nchmx_loc),rt(nchmx_loc),w(*)
      complex s(*),tfunc(3,nchmx_loc),x(nchmx_loc),h(nchmx_loc),
     & xx(nchin,*)
 
      if(lref) then
         nch = nchstk - 2
      else
         nch = nchstk
      endif

      if(nf.lt.2) then
         print*,'not enough data to process'
         return
      end if

      nch2 = nch-2
      n2 = nch/2
      ns = (nchstk*(nchstk+1))/2

c........ Weighting by power; not used much (if at all)
ccc      2-26-98   Removed; Egbert
c      if(wtpram.gt.1.e-20) then
c         print*,'weighting by power'
c         call wtbpwr(xx,nchin,nf,wtpram)
c      end if

c>>>>>>>>>>this block is executed for non-robust Tfs<<<<<<<<<<<<
      if((.not.lrbst) .or. (nf.lt.3)) then
c.........zero and then stack (weighted) sdm
         do i=1,ns
            s(i)=(0.,0.)
         enddo
         do i = 1,nf
            call stack(xx(1,i),nchstk,s)
         enddo
         nu = nf
         do i=1,nch2
            ru(i)=nf*2.
         enddo
         return
      end if
c>>>>>>>>>>>>> end of non-robust block<<<<<<<<<<<<<<<<<<<<<<<<<<<

c
c>>>>>>>>>>>> ROBUST TF BLOCK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c                      (includes rest of subroutine)

c.......estimate horizontal mag field sdm robustly
      call mkrhat(edfl1,nf,xx,nchin,h,ierr,w)
      if(ierr .eq. -1) then
         print*,'error in mkrhat'
         return
      end if

c.........zero sdm, set up leverage point parameters
      do i=1,ns
         s(i)=(0.,0.)
      enddo
      p1 = float(nf)**alpha
      p2 = c2*p1
      p1 = c1*p1

c.....compute weights for reducing effect of leverage points
      call edfwts(xx,nchin,nf,p1,p2,h,w)

c>>>>>>>>>>>>>>>> initial tf estimate (weighted ls)
      iter=0

c.....multiply data by weights and form sdm for nchstk channels
      nu = 0
      do i = 1,nf
         if(w(i) .gt. 0.)then
            nu = nu + 1
            if(w(i) .lt. 1.) then
               do j= 1,nchstk
                   xx(j,i) = w(i)*xx(j,i)
               enddo
            end if
            call stack(xx(1,i),nchstk,s)
         end if
      enddo

c.....compute initial tf estimate
cnew      print*,'initial estimate,;nf,nu:',nf,nu
      if(nu .lt. 2) return
      if(lref) then
         call trlsrr(s,nch,nu)
      else
         call tranls(s,nch,nu)
      end if

c.........save initial tf estimates, error scales
      call savtf(s,nch,tfunc)
      do i=1,nch2
          ldone(i)=.false.
         scale(i)=sqrt(real(tfunc(3,i))/2.)
      enddo
 
c>>>>>>>> Robust tf estimates:
c                loop with non redescending influence curve
20    continue
      iter=iter+1
c.....zero sdm
      do i = 1,ns
         s(i) = 0
      enddo
      nu = 0
      do ifreq=1,nf
         if(w(ifreq) .gt. 0.) then
c.......... data point is possibly useful; move into temporary array
            do i = 1,nchstk
               x(i) = xx(i,ifreq)
            enddo
                                        
c..........modify data if necessary;(pull outliers toward predicted values)
            do i = 3,nch
               i2 = i-2
               call modobs(x,x(i),wt,tfunc(1,i2),scale(i2),rus(i2),
     &              rt(i2))
            enddo
            
c.......... compute median residual power for nch2 channels
            call sort(rt,nch2)
            if(rt(n2).lt.c3) then
c..............stack modified data
               nu = nu + 1
               call stack(x,nchstk,s)
            end if
         end if
      enddo    ! ifreq
 
c.....compute new transfer function estimates
      if(nu .lt. 2) then
         print*,'threw all data away!!!  iteration = ',iter,'  nf = '
     &         ,nf,'nu = ',nu
         return
      end if
      if(lref) then
         call trlsrr(s,nch,nu)
      else
         call tranls(s,nch,nu)
      end if 

c.....test for convergence of TFs for each predicted channel
      do i=3,nch
         is=sindex(i,1)
         tst=(abs(s(is)-tfunc(1,i-2)))**2
         is=is+1
         tst=tst+(abs(s(is)-tfunc(2,i-2)))**2
         tst=sqrt(tst/(abs(tfunc(1,i-2))**2+abs(tfunc(2,i-2)**2)))
         ldone(i-2)=(tst.le.tol)
      enddo
 
c.....save current tf estimates, error scales
      call savtf(s,nch,tfunc)
      do i=1,nch-2
         scale(i)=sqrt(real(tfunc(3,i))/(varcor*2.))
      enddo
 
c.....if not converged for all channels go to top of loop
      lastit=.true.
         do i=1,nch-2
            lastit=lastit.and.ldone(i)
         enddo
      if((.not.lastit).and.(iter.lt.itmax)) go to 20
      if(iter.eq.itmax) print*,'maximum number of iterations reached'
        
c>>>>>>>>> if all channels have converged (or itmax reached):
c.....do final iteration(s) using redescending influence curve
 
c.....set up for final iteration
      nu = 0
      do i=1,ns
         s(i)=(0.,0.)
      enddo
      do i = 1,nch-2
         ru(i) = 0.
      enddo
      do ifreq=1,nf
         if(w(ifreq).gt.0.) then
c.......... move record into temporary storage
            do i = 1,nchstk
               x(i) = xx(i,ifreq)
            enddo
c...........modify data; pull outliers toward predicted values
            do i = 3,nch
               i2 = i-2
               call modobs(x,x(i),rdcndwt,tfunc(1,i2),scale(i2),rus(i2)
     &            ,rt(i2))
            enddo
c.......... compute median residual power (i.e. median for 3 channels)
            call sort(rt,nch2)
            if(rt(n2).lt.c3) then
c.....................stack modified data
               nu = nu + 1
               call stack(x,nchstk,s)
               do j=1,nch2
                  ru(j) = ru(j) + rus(j)
               enddo
            end if
         end if
      enddo

cnew      print*,'final result; iter,nf,nu ',iter,nf,nu
      return
      end
