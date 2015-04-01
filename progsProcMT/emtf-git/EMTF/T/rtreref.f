c____________________________________________________________
c 
      subroutine rtreref(z,nbt,nz,ldf,nch,nsmx)
c
ccc   computes addmitance estimates and error covariance
ccc   then converts to impedance estimate with errors
ccc   Overwrites results
c     in array z; output can then be used in MT interpretaion
c          parameter routines exactly as usual after call to tranls

c      NOTE: does not currently handle Z transfer functions
c      properly ... these are computed in a different call to gdsint

      complex y(2,2),eei(2,2),rsig(2,2),z(nsmx,*),a(2,2),b(2,2),det
     & ,e(2,2),h(2,2)
      real coh(2)
      integer ldf(*)               
 
      if(nch.eq.5) then
         ih1 = 10
         ih2 = 14
         ih3 = 15
         ix1 = 7
         ix2 = 8
         ix3 = 11
         ix4 = 12
      else if(nch.eq.4) then
         ih1 = 6
         ih2 = 9
         ih3 = 10
         ix1 = 4
         ix2 = 5
         ix3 = 7
         ix4 = 8
      endif
      do ib = 1,nbt
ccc      save <ee>
         e(1,1) = z(1,ib)
         e(1,2) = z(2,ib)
         e(2,1) = conjg(z(2,ib))
         e(2,2) = z(3,ib)
ccc      and <hh>
         h(1,1) = z(ih1,ib)
         h(1,2) = z(ih2,ib)
         h(2,1) = conjg(z(ih2,ib))
         h(2,2) = z(ih3,ib)
ccc      compute admittance     
         call tranls(z(1,ib),nch,ldf(ib))
         z(ih2,ib) = conjg(z(ih2,ib)) - z(ix3,ib)*conjg(z(ix1,ib))
     &        *e(1,1) - z(ix4,ib)*conjg(z(ix1,ib))*e(2,1)
     &        - z(ix3,ib)*conjg(z(ix2,ib))*e(1,2)
     &         - z(ix4,ib)*conjg(z(ix2,ib))*e(2,2)
      
ccc      invert estimated admittance tensor
         det = z(ix1,ib)*z(ix4,ib)-z(ix2,ib)*z(ix3,ib)
         y(1,1) = z(ix4,ib)/det
         y(2,2) = z(ix1,ib)/det
         y(1,2) = -z(ix2,ib)/det
         y(2,1) = -z(ix3,ib)/det

         eei(1,1) = z(1,ib)
         eei(1,2) = conjg(z(2,ib))
         eei(2,1) = z(2,ib)
         eei(2,2) = z(3,ib)
         rsig(1,1) = real(z(ih1,ib))*(ldf(ib)-2)
         rsig(1,2) = conjg(z(ih2,ib))
         rsig(2,1) = z(ih2,ib)
         rsig(2,2) = real(z(ih3,ib))*(ldf(ib)-2)

ccc      compute coherence
         do i = 1,2
            coh(i) = 0.
            do k = 1,2
               do l = 1,2
                  coh(i) = coh(i) + y(i,k)*conjg(y(i,l))*h(k,l)
               enddo
            enddo
            coh(i) = e(i,i)/coh(i)
            coh(i) = min(coh(i),.999)
            coh(i) = max(0.,coh(i))
         enddo

         do i = 1,2
            do j = 1,2
               a(i,j) = (0.,0.)
               b(i,j) = (0.,0.)
               do k = 1,2
                  do l = 1,2
                     a(i,j) = a(i,j) + y(k,i)*conjg(y(l,j))*eei(k,l)
                     b(i,j) = b(i,j) + y(i,k)*conjg(y(j,l))*rsig(k,l)
                  enddo
               enddo
            enddo
         enddo

         z(1,ib) = a(1,1)
         z(3,ib) = a(2,2)
         z(2,ib) = a(1,2)
         z(ih1,ib) = cmplx(real(b(1,1))/(ldf(ib)-2),real(coh(1)))
         z(ih3,ib) = cmplx(real(b(2,2))/(ldf(ib)-2),real(coh(2)))
         z(ih2,ib) = b(1,2)/(ldf(ib)-2)
         z(ix1,ib) = y(1,1)
         z(ix2,ib) = y(1,2)
         z(ix3,ib) = y(2,1)
         z(ix4,ib) = y(2,2)
      enddo    ! ib
      return
      end
