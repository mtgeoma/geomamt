      subroutine rotsdm(zz, u, nch)
c     rotsdm rotates the spectral density matrix applying a rotation matrix u
c     zz is the spectral density matrix in compact storage mode, i. e. 
c     only the diagonal terms and those below the diagonal are stored. There are
c     nch*(nch+1)/2 elements.
c     u is the full rotation matrix, i.e a nch x nch matrix, though many elements
c     are 0.

      complex zz(*), zzz(7,7), a(7,7)
      real u(7,7)


c     make the full matrix zzz from the compact SDM zz
      kj = 0
      do k = 1, nch
         do j = 1, k
            kj = kj + 1
            zzz(j,k) = (zz(kj))
            zzz(k,j) = conjg(zz(kj))
         enddo
      enddo
         
c     now multiply u from left with zzz
      do i = 1, nch
         do j = 1,nch
            a(i,j) = (0.0,0.0)
            do k = 1,nch
               a(i,j) = a(i,j) + u(i,k)*zzz(k,j)
            enddo
         enddo
      enddo
c     and now the transposed of u from the right to a
      do i = 1, nch
         do j = 1,nch
            zzz(i,j) = (0.0,0.0)
            do k = 1,nch
               zzz(i,j) = zzz(i,j) + a(i,k) * u(j,k)
            enddo
         enddo
      enddo

c     finally, restore the full matrix zzz in the compact zz
      kj = 0
      do k = 1, nch
         do j = 1,k
            kj = kj + 1
            zz(kj) = zzz(k,j)
         enddo
      enddo

      return
      end

