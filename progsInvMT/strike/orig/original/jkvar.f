      subroutine jkvar( x, n, var )
      implicit none
      integer n
      real x(n), var, xn

      real xtot, mean, d1mean
      integer i

c...derive jackknife estimate of variance of mean of data series X(N)
c
c   jknmean = n*mean - (n-1)/n * sum[ (delete-1 means) ]
c             (for a linear system, this is a straight mean)
c
c   jknvar = (n-1)/n * sum[ ( (delete-1 means) - mean )**2 ]

c...check n > 1
      if( n.le.1 ) then
        write(*,'(a,i5)')'jkvar: n too small: n = ', n
        var = -1.
        return
      endif

      xn = float(n)

c...get total sum
      xtot = 0.0
      do i = 1, n
        xtot = xtot + x(i)
      enddo
      mean = xtot / xn

c...sum delete-1 estimates of variance
      var = 0.0
      do i = 1, n
        d1mean = (xtot - x(i))/(xn-1.)

        var = var + (d1mean - mean)*(d1mean - mean)
      enddo


      var = ( (xn-1.)/xn ) * var

      return
      end
