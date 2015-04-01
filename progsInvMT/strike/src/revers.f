      subroutine revers( n, a )

      implicit none

      real a(*), dummy
      integer i, n

c     subroutine to reverse real array "A" of length "N"

      do i = 1, n/2
        dummy    = a(i)
        a(i)     = a(n-i+1)
        a(n-i+1) = dummy
      enddo

      return
      end
