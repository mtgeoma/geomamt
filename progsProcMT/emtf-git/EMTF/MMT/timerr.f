      subroutine timerr(nsta,nt,ntape,nch,nd,orient,ih,ie,twin,xx,yy,
     &    lx,isuse,iterrb,nterrb,terr,sig,u,wrk,cwrk,iwrk,s)

c      Computes timing error corrections :  New version 12-10-93
c        IF twin is length of window in seconds,
c       timerr is time shift DT in seconds

      include 'iosize.inc'

      parameter (pi = 3.14159, pi2 = pi/2.)

      integer iterrb(3,nterrb),isuse(2,*),nch(nsta),ntape(nsta),ih(*),
     & iband(2),iwrk(*),ie(*),orient(2,nt)
      real terr(nsta),tx(1000),ev(100),ty(1000),txx(1000),sig(*),
     & cov(100),wt(1000),declt(100),wrk(*),b(500),sdiag(100),theta
      complex c(2),a(2),xx(*),s(*),cwrk(*),u(nt,*),y(500)
      logical lquad,lx(*),llx,lch(500),e_and_h

      do 5 i = 1,nsta
5     terr(i) = 0.
      idec = 1

      id = 1
      ii = 0
      nnd = 0

      write(*,*) 'Timing Error Bands', iterrb(1,1), iterrb(2,1)
     &    ,iterrb(3,1)
      write(*,*) 'ih',(ih(k),k=1,3)
      do 100 ib = 1,nterrb
         do 80 i = iterrb(1,ib),iterrb(2,ib),iterrb(3,ib)
         iband(1) = i
         iband(2) = i + iterrb(3,ib) - 1
         ii = ii + 1

ccc         make array data records ( ixs,ixf are dummies here)
ccc         get only data with all stations
         llx = .false.
         call mkrec(nd,id,iband,nsta,ntape,nch,nt,isuse,
     &  nfreq,xx,llx,lx,ixs,ixf,cwrk,iwrk)
ccc    to use all data in stacking routine need to set lx = .ture.
ccc    for all data points
         do ifreq=1,nfreq
            lx(ifreq) = .true.
         enddo
         do ich = 1,nt
            lch(ich) = .true.
         enddo 

         if(nfreq .le. 2) then
            wt(ii) = 0.
            go to 80
         end if

ccc      calling rbstk with itmax=0 makes this a standard unweighted sdm           
         write(*,*) 'nt,nfreq',nt,nfreq
         iter_max = 10
         call rbstk(xx,nt,nt,nfreq,nfreq,s,yy,iter_max,lx,lch)
         icon = 1
         call diag(s,nt,sdiag) 
         call rsp(s,nt,2,u,ev,icon,wrk,sdiag)
         write(*,*) 'ev  = ',(ev(k),k=1,10)
ccc      change coordinates to a common geographic system
         e_and_h = .false.
         theta = 0.
         call cc_geog(u,nt,2,nsta,ih,ie,orient,theta,e_and_h)


ccc      use total power to detrmine weights for LS fit of timing error parameters
         period = ii
         wt(ii) = 0.
            do 45 j = 3,nt
            wt(ii) = wt(ii) + ev(j)
45          continue
            wt(ii) = wt(ii)/((nfreq-2)*(ev(1)+ev(2)))
            wt(ii) = 1./sqrt(wt(ii))
ccc            wt(ii) = 1.
          
ccc        find response vector for unit H at site 1
         call mkbih(b,nt,ih,ie,nsta,1)
         c(1) = (1.,0.)
         c(2) = (0.,0.)
         call anfld(u,nt,nt,b,c,a,y)
           print*,y(1),y(6)
         
         do j = 2,nsta
            nnd = nnd+1
            tx(nnd) = atan( aimag(y(ih(j))) / real( y(ih(j)) )  )

            if( real( y(ih(j)) ) .lt. 0.) then
               if(aimag(y(ih(j))) .lt. 0.) then
                  tx(nnd) = tx(nnd) - pi
               else
                  tx(nnd) = tx(nnd) + pi
               end if
            end if
         enddo

         c(1) = (0.,0.)
         c(2) = (1.,0.)
         call anfld(u,nt,nt,b,c,a,y)
           print*,y(2),y(7)

         do j = 2,nsta
            nnd = nnd+1
            tx(nnd) = atan( aimag(y(ih(j)+1)) / real( y(ih(j)+1) )  )
            if( real( y(ih(j)+1) ) .lt. 0.) then
               if(aimag(y(ih(j)+1)) .lt. 0.) then
                  tx(nnd) = tx(nnd) - pi
               else
                  tx(nnd) = tx(nnd) + pi
               end if
            end if
         enddo

80       continue
100   continue


      print*, 'in terr ; nsta = ',nsta
      np = 2*nterrb+1
         do 200 ista = 2,nsta

c      this block takes care of the problem that might arise if
c       atan chooses the wrong quadrant
c         Basic idea: first difference succcesive phases; if
c          phase changes by more than pi/2 assume that the wrong choice
c       been made, and hence take the other branch;  should work fine
c     unless timing errors exceed several times the sampling frequency.

         nnd2 = nnd/2
         do i = 1,nnd2
            ty(i) = tx(2*(nsta-1)*(i-1)+ista-1)
         enddo

         lquad = .false.
         do 117 i = 1,nnd2-1
            tyi = ty(i+1) - ty(i)
            if(abs(tyi).gt.pi)then
               lquad = .true.
               if(tyi.gt. 0.) then
                  ty(i) = tyi - 2.*pi
               else
                  ty(i) = 2.*pi + tyi
               end if
            else
               ty(i) = tyi
            end if
117      continue

         if(lquad) then
            do 118 i = 2,nnd2
            i0 = 2*(nsta-1)*(i-2)+ista-1
            i1 = i0 + 2*(nsta-1)
            tx(i1) = tx(i0)+ty(i-1)
118         continue
         end if

            do 216 i = 1,nnd2
216         ty(i) = tx((2*i-1)*(nsta-1)+ista-1)
         lquad = .false.
            do 217 i = 1,nnd2-1
            tyi = ty(i+1) - ty(i)
            if(abs(tyi).gt.pi)then
               lquad = .true.
               if(tyi.gt. 0.) then
                  ty(i) = tyi - 2.*pi
               else
                  ty(i) = 2.*pi + tyi
               end if
            else
               ty(i) = tyi
            end if
217         continue
         if(lquad) then
            do 218 i = 2,nnd2
            i0 = (2*i-3)*(nsta-1)+ista-1
            i1 = i0 + 2*(nsta-1)
            tx(i1) = tx(i0)+ty(i-1)
218         continue
         end if

c            print*,'station',ista
c            print*,'H phases (unweigthed)'
c            print*,(2*i+20,tx(2*(nsta-1)*(i-1)+ista-1),i=1,6)
c            print*,'D phases (unweigthed)'
c            print*,(2*i+20,tx(2*(nsta-1)*(i-1)+nsta+ista-2),i=1,6)

         nnd = 0
         jj = ista-1
         ii = 0
            do 120 ib = 1,nterrb
               dt = float(iterrb(3,ib))/2.
               do 120 i =  iterrb(1,ib),iterrb(2,ib),iterrb(3,ib)
               ii = ii + 1
               if(wt(ii) .gt. 1.0e-10) then
                  nnd = nnd+1
                  ty(nnd) = wt(ii) *tx(jj)
                  txx(nnd) = 2.*pi*wt(ii)*(i+dt)/twin
                  nnd = nnd+1
                  jj = jj+nsta-1
                  ty(nnd) = wt(ii)*tx(jj)
                  txx(nnd) = 2.*pi*wt(ii)*(i+dt)/twin
                  jj = jj+nsta-1
               end if
120         continue

         ll = nnd
         do 150 ib = 1,nterrb
         ii = 0
            do 130 jb = 1,nterrb
            do 130 j = iterrb(1,jb),iterrb(2,jb),iterrb(3,ib)
            ii = ii + 1
            if(wt(ii) .gt. 1.0e-10) then
               if(ib.eq.jb) then
                  ll = ll + 1
                  txx(ll) = wt(ii)
                  txx(ll+nnd) = 0.
                  ll = ll + 1
                  txx(ll) = 0.
                  txx(ll+nnd)= wt(ii)
               else
                  ll = ll + 1
                  txx(ll) = 0.
                  txx(ll+nnd) = 0.
                  ll = ll + 1
                  txx(ll) = 0.
                  txx(ll+nnd) = 0.
               end if
            end if

130         continue
            ll = ll + nnd

150      continue

         call lsv(txx,ty,nnd,np,terr(ista),cov,wrk)
         sig(ista) = sqrt(cov(1))
           
         print*,'terr',terr(ista)
         print*,'sig',sig(ista)

c         if(abs(terr(ista)).lt.(2.*sig)) terr(ista) = 0.
200      continue
         return
         end
c______________________________________________________________________
c
      subroutine terrfix(xx,nt,nf,nsta,ih,period,terr)
      complex xx(nt,nf),terrc(100)
      real terr(nsta)
      integer ih(*)
      parameter (pi=3.14159)

c       print*,'nt,nf',nt,nf
c       print*,'terr',terr
c       print*,'period',period
c       print*,'nsta',nsta
c       print*,'ih',(ih(ista),ista=1,nsta+1)
      do ista=1,nsta
         nch = ih(ista+1)-ih(ista)
         t = 2.*pi*terr(ista)/period
         do k=1,nch
            i = ih(ista)+k-1
            terrc(i)=cmplx(cos(t),-sin(t))
         enddo
      enddo

      do i = 1,nf
         do j = 1,nt
            xx(j,i) = xx(j,i)*terrc(j)
         enddo
      enddo
      return
      end
