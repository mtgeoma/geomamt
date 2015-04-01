c______________________________________________________________________
c
      subroutine filtpc(nt,nev,u,ev,u2,ev2,wrk,ih,nsta)

      include 'nstamx.inc'

      complex u(nt,nev),u2(nt,nev),wrk(*)
      real ev(nev),ev2(nev),eh(ntmx),hh,ee
      integer ih(nsta)

ccc      given principle components of a scaled SDM
ccc          form the "truncated SDM", of rank nev
ccc     (scaled back into original coordinate system)
ccc     (i.e., truncate to the dominant nev PCs and reconstruct the SDM)
ccc      and then redo the SDM to find orthogonal evecs, and meaningful
ccc      estimates of signal power.

cc     For E fields, we want to use the impedance to put everything in
ccc     the same physical units ... will use a crude estimate of impedance first

ccc      reconstruct truncated SDM
      call sing_s(nt,nev,ev,u,wrk)
ccc      crude estimate of array/component averaged impedance amplitude
      hh = 0.0
      ee = 0.0
      do ista = 1,nsta
         ii = (ih(ista)*(ih(ista)+1))/2
         hh = hh + wrk(ii)
         ii = ii + ih(ista)+1
         hh = hh + wrk(ii)
         ii = (ih(ista+1)*(ih(ista+1)-1))/2
         ee = ee + wrk(ii)
         ii = ii - ih(ista+1) + 1
         ee = ee + wrk(ii)
      enddo
      do i = 1,nt
         eh(i) = 1.0
      enddo
      do ista = 1,nsta
         eh(ih(ista+1)-2) = sqrt(hh/ee)
         eh(ih(ista+1)-1) = sqrt(hh/ee)
      enddo
      ij = 0
      do i = 1,nt
         do j = 1,i
            ij = ij + 1
            wrk(ij) = wrk(ij)*eh(i)*eh(j)
         enddo
      enddo

ccc     redo SDM
      icon = 0
      nwrk1 = (nt*(nt+1))/2 + 1
      call rsp(wrk,nt,nev,u2,ev2,icon,wrk(nwrk1))

ccc      change E fields back to mV/Km
      do i = 1,nt
         do j = 1,nev
            u2(i,j) = u2(i,j)/eh(i)
         enddo
      enddo

ccc      change phase to minimize phase of horizontal magnetics
      call chngph(u2,nt,ih,nsta,nev)

      return
      end
c______________________________________________________________________
c
      subroutine chngph(b,nt,ih,nsta,nvec)
c     change phases of multiple station vector to make
c       rms imaginary part of horizonatal magnetics zero
      complex tc,ac
      real b(2,nt,nvec)
      integer ih(*)

      do l = 1,nvec
         sxx = 0.
         sxy = 0.
         syy = 0.
         sy = 0.
            do i = 1,nsta
            sy = sy + b(2,ih(i),l)
            do k = 0,1
               sxx = sxx + b(1,ih(i)+k,l)*b(1,ih(i)+k,l)
               sxy = sxy + b(1,ih(i)+k,l)*b(2,ih(i)+k,l)
               syy = syy + b(2,ih(i)+k,l)*b(2,ih(i)+k,l)
            enddo
         enddo

         r = 2.*sxy/(sxx-syy-sqrt((sxx-syy)**2+4.*sxy*sxy))
         bb = 1./sqrt(1+r*r)
         aa = r*bb
         if(sy.gt.0) then
            bb = -bb
            aa = -aa
         end if
         ac = cmplx(aa,bb)

         do i = 1,nt
            tc = cmplx(b(1,i,l),b(2,i,l))
            tc = tc*ac
            b(1,i,l) = real(tc)
            b(2,i,l) = aimag(tc)
         enddo
      enddo   ! l=1,nvec

      return
      end
c______________________________________________________________________
c
      subroutine sing_s(nt,nev,ev,u,ss)
      complex ss(*),u(nt,nev)
      real ev(nev)

      ij = 0
      do i = 1,nt
         do j = 1,i
            ij = ij + 1
            ss(ij) = (0,.0)
            do k = 1,nev
               ss(ij) = ss(ij) 
     &                  + ev(k)*u(i,k)*conjg(u(j,k))
            enddo
         enddo
      enddo
      return
      end
