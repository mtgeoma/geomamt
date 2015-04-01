c**********************************************************
c   CHI SQUARE ERROR 
c
c     calculate the chi-squared error of misfit between
c     the data (zdata) and the response (zresp), normalized by
c     the variance (squared standard error) of the data (std_dev) 
c
c***********************************************************      

      subroutine calc_error(zdata,zresp,std_dev, error)

      implicit none

      complex*16 zdata(2,2),zresp(2,2)
      real*8     error, delta
      real*8     std_dev(2,2)
      integer    i, j

      error = 0.0D+00
      do i = 1, 2
        do j = 1, 2
          delta = cdabs( zdata(i,j)-zresp(i,j) )
          error = error + ( delta**2 )/( std_dev(i,j)**2 )
        enddo
      enddo

      return
      end
