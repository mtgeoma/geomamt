c
c***********************
c
      subroutine autocor(x,npts,nlag,nch,ich,r,rr)                
c
c       computes autocorrlation matrix for a segment
c      of time series x; auto-correlations up to lags
c       of nlag are returned in r;
c       a guarranteed to be pos. semi-def. (but non-Toeplitz)
c       covariance matrix out to lag nlag-1 is returned in rr

      real x(nch,npts),r(*),rr(*)
c         damp controls damping of covariance matrix to increase
c        stability of pre-whitening and avoid the introduction of
c        spurious peaks in the whitened spectra

      parameter (damp = .01)

      do 20 j=1,nlag+1
      r(j)=0.         
         do 10 k=1,npts
         r(j)=r(j)+x(ich,k)*x(ich,k+j-1)
10       continue
20    continue
 
c      make Toeplitz matrix from autocovariance
      ii = 0
      do 30 i=1,nlag
      do 30 j = 1,i
      ii = ii + 1
      rr(ii) = r(i-j+1)
30    continue

c      modify Toeplitz matrix to make pos. def matrix
      do 40 j=2,nlag
c       modify jth column
         do 40 i=j,nlag
         ij = (i*(i-1))/2+j
            do 35 k = 1,j-1
            rr(ij) = rr(ij) - x(ich,k)*x(ich,k+i-j)
     &                            + x(ich,npts+k)*x(ich,npts+k+i-j)
35          continue
40    continue

      do 45 j = 1,(nlag*(nlag+1))/2
45    rr(j) = rr(j)/npts

      do 50 j =1,nlag+1
50    r(j) = r(j)/npts
 
c       damp autocov coefficient matrix

      do 60 i = 1,nlag
      ii = (i*(i+1))/2
      rr(ii) = rr(ii) + r(1)*damp
60    continue
      return
      end 
