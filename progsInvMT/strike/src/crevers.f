      subroutine crevers( n, a )

      implicit none

      complex a(*), dummy
      integer i, n

c     subroutine to reverse complex array "A" of length "N"

      do i = 1, n/2
        dummy    = a(i)
        a(i)     = a(n-i+1)
        a(n-i+1) = dummy
      enddo

      return
      end
