c***********************************************************************
c
c  subroutine extreme - determines lower and upper error bar
c
c***********************************************************************

      subroutine extreme(vect,len, vmax, vmin )

      real vect(len), vmax, vmin, x
      integer len

      vmax=vect(1)
      vmin=vect(1)
            
      do i = 2, len
          x = vect(i)
          if( x. gt. vmax) vmax = x
          if( x .lt. vmin) vmin = x
      enddo

      return
      end
