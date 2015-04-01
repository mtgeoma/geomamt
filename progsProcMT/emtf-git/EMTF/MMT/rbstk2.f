      subroutine rbstk2(xx,nt,nf,s,y,itmax)

ccc   new version assumes all data is present, with each channel a
ccc   a column of the data matrix XX; no logical masks here
      include 'nstamx.inc'
      parameter (tol=.001,epsilon = 1.e-4,r0fac=1.5)
      complex xx(nf,nt),s(*),y(nt,nf),y1(ntmx),cdot

      logical lrbst

      lrbst = (itmax.gt.0) .and. (nf .gt. 2*nt)

ccc   zero sdm array S
      r0 = r0fac*float(nt)
      ns = ((nt+1)*nt)/2
      do ij = 1,ns
         s(ij) = (0.,0.)
      enddo

ccc   copy data to use into array y and compute initial SDM
ccc   note: different row-column convention for Y and X
      w = 1.0
      ii = 0
      do i = 1,nf
         do j=1,nt
            y(j,i) = xx(i,j)
         enddo
         call wstack(y(1,i),w,nt,s)
      enddo

      do ij = 1,ns
         s(ij) = s(ij)/nf
      enddo

cccc  Raw stack is done, and Y contains a copy of data matrix
cccc  if itmax = 0, don't do robust iterations
      if(.not.lrbst) return

ccc     iterative loop to calculate robust stacking weights
      do 100 iter = 1,itmax
c            form cholesky decomp of current "normalized SDM"
         call stabls(s,nt,epsilon)
         call cchdeci(s,nt,ier)
         if(ier.eq.-1) then
             write(*,*) '!!!!!!! ERROR IN RBSTK : S SINGULAR!!!!!'
         end if
c             normalize data vectors
         call cltslv(s,nt,y,nf)
c             compute weights and stack
         wtot = 0.
         do 40 ij = 1,ns
40          s(ij) = (0.,0.)
         do 50 i=1,nf
            r = cdot(y(1,i),y(1,i),nt)
            if(r.gt.r0) then
               w = r0/r
c               print*,'r,w',r,w
            else
               w = 1.0
            end if
            wtot = wtot + w
            call wstack(y(1,i),w,nt,s)
50          continue
c            print*,'wtot = ',wtot

         schk = 0.0
         ij = 0
         do 60 i = 1,nt
         do 60 j = 1,i
         ij = ij + 1
            s(ij) = s(ij)/wtot
c            print*,'i,j,ij,s(i,j)',i,j,ij,s(ij)
            if(i.eq.j) then
               schk = schk + (abs(1.-s(ij)))**2.
            else
               schk = schk + (abs(s(ij)))**2.
            end if
60          continue
c    check to see if normalized weigthed sdm is near identity
         if(schk.le.tol) go to 200 
100      continue 
         write(*,*) 'Maximum number of iterations exceeded'
         write(*,*) 'itmax,schk,tol = ',itmax,schk,tol

ccc        end of iterative computation of weights
200   continue
c        
      do ij = 1,ns
         s(ij) = (0.,0.)
      enddo

      ii = 0
      wtot = 0.
      do i=1,nf
         r = cdot(y(1,i),y(1,i),nt)
         if(r.gt.r0) then
            w = r0/r
         else
            w = 1.0
         end if
         wtot = wtot + w
         do j = 1,nt
            y1(j) = xx(i,j)
         enddo
         call wstack(y1,w,nt,s)
      enddo

      do ij = 1,ns
         s(ij) = s(ij)/wtot
      enddo

      return
      end
c______________________________________________________________________
c
        subroutine wstack(x,w,nch,s)
        complex x(nch),s(*),y
        real w
        ij=0
        do 10 i=1,nch
        y=w*x(i)
        do 10 j=1,i
           ij=ij+1
10         s(ij)=s(ij)+y*conjg(x(j))
        return
        end
c______________________________________________________________________
      subroutine stabls(s,nt,epsilon)

      complex s(*)
      real trs,epsilon

c      compute trace ....
      trs = 0.0
      ij = 0
      do 35 i = 1,nt
         ij = ij + i
         trs = trs + s(ij)
35       continue

      trs = epsilon*trs/nt
      ij = 0
      do 36 i = 1,nt
         ij = ij + i
         s(ij) = s(ij) + trs
36       continue
      return
      end
c______________________________________________________________________
c
      subroutine rbstk(xx,nt_tot,nt,nf_tot,nf,s,y,itmax,luse,lch)

      include 'nstamx.inc'
      parameter (tol=.001,epsilon = 1.e-4,r0fac=1.5)
      complex xx(nt_tot,nf_tot),s(*),y(nt,nf),y1(ntmx),cdot

ccc     luse is a logical mask which allows for stacking of subsets of data
ccc     lch is a logical mask which allows for eliminating components
ccc        ( or whole stations) from the SDM
ccc        this mask should have exactly nt .true. elements
ccc      Use these masks for doing sub-arrays : in this case array xx will contain
ccc       all data; luse and lch will tell which segments/stations to use

      logical luse(nf_tot),lch(nt_tot),lrbst

ccc      first check to see if lch and nt are compatable

      lrbst = (itmax.gt.0) .and. (nf .gt. 2*nt)
      nt0 = 0
      do i = 1,nt_tot
         if(lch(i)) nt0 = nt0+1
      enddo
      if(nt0.ne.nt) then
         write(*,*) 'Mask lch and number of components nt',
     &    ' in rbstk are incompatable; STOPPING'
         stop
      endif

ccc      zero sdm array S
      r0 = r0fac*float(nt)
      ns = ((nt+1)*nt)/2
      do ij = 1,ns
         s(ij) = (0.,0.)
      enddo

ccc     copy data to use into array y and compute initial SDM
      w = 1.0
      ii = 0
      do i = 1,nf_tot
         if(luse(i)) then
            ii = ii + 1
            jj = 0
            do j=1,nt_tot
               if(lch(j)) then
                  jj = jj + 1
                  y(jj,ii) = xx(j,i)
               end if
            enddo
            call wstack(y(1,ii),w,nt,s)
         endif
      enddo
      nf = ii
      do ij = 1,ns
         s(ij) = s(ij)/nf
      enddo

cccc    Raw stack is done, and Y contains data to use for array
ccc      nf is number of frequencies to use
ccc     NOTE: during iterations to follow, the original data in Y is lost;
ccc        still need to use data from XX to compute SDM estimate

cccc     if itmax = 0, don't do robust iterations
      if(.not.lrbst) return

ccc     iterative loop to calculate stacking weights
      do 100 iter = 1,itmax
c            form cholesky decomp of current "normalized SDM"
         call stabls(s,nt,epsilon)
         call cchdeci(s,nt,ier)
         if(ier.eq.-1) then
             write(*,*) '!!!!!!! ERROR IN RBSTK : S SINGULAR!!!!!'
         end if
c             normalize data vectors
         call cltslv(s,nt,y,nf)
c             compute weights and stack
         wtot = 0.
         do 40 ij = 1,ns
40          s(ij) = (0.,0.)
         do 50 i=1,nf
            r = cdot(y(1,i),y(1,i),nt)
            if(r.gt.r0) then
               w = r0/r
c               print*,'r,w',r,w
            else
               w = 1.0
            end if
            wtot = wtot + w
            call wstack(y(1,i),w,nt,s)
50          continue
c            print*,'wtot = ',wtot

         schk = 0.0
         ij = 0
         do 60 i = 1,nt
         do 60 j = 1,i
         ij = ij + 1
            s(ij) = s(ij)/wtot
c            print*,'i,j,ij,s(i,j)',i,j,ij,s(ij)
            if(i.eq.j) then
               schk = schk + (abs(1.-s(ij)))**2.
            else
               schk = schk + (abs(s(ij)))**2.
            end if
60          continue
c    check to see if normalized weigthed sdm is near identity
         if(schk.le.tol) go to 200 
100      continue 
         write(*,*) 'Maximum number of iterations exceeded'
         write(*,*) 'schk,tol = ',schk,tol

ccc        end of iterative computation of weights
200   continue
c        
      do ij = 1,ns
         s(ij) = (0.,0.)
      enddo

      ii = 0
      wtot = 0.
      do i=1,nf_tot
         if(luse(i)) then
            ii = ii + 1
            r = cdot(y(1,ii),y(1,ii),nt)
            if(r.gt.r0) then
               w = r0/r
            else
               w = 1.0
            end if
            wtot = wtot + w
            jj = 0
            do j = 1,nt_tot
               if(lch(j)) then 
                  jj = jj + 1
                  y1(jj) = xx(j,i)
               endif
            enddo
            call wstack(y1,w,nt,s)
         endif
      enddo

      do ij = 1,ns
         s(ij) = s(ij)/wtot
      enddo

      return
      end
