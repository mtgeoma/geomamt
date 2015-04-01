c
c**************************************
c
      subroutine cldrft(x,nfreq,nch,dr,nwin,cdb,tset)
     
c    subroutine corrects for linear clock drift
c     clock time is tc; true time is t
c     then assuming tc = t + cdb*t  this routine corrects for
c     the phase shift in the complex fourier coefficients for
c     a window centered at time tset (where tset is in seconds after
c    clock reset) it is assumed that at t=0 (i.e. at clock reset)
c    that the clock time is correct
       
      complex tc,t
      real x(nch,2,nfreq)
      parameter(pi2 = 6.28318)                      

      t1 =  - pi2*cdb*tset/(dr*nwin)

      do 10 i = 1,nfreq
      tc = cmplx(0.,i*t1)
      tc = cexp(tc)
         do 5 j = 1,nch
         t = tc*cmplx(x(j,1,i),x(j,2,i))
         x(j,1,i) = real(t)
5        x(j,2,i) = aimag(t)
10    continue

      return
      end
