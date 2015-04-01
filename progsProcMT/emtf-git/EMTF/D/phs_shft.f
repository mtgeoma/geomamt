ccc_____________________________________________________________________
ccc
      subroutine phs_shft(dt,x,nfreq,nch,dr,nwin)
     
ccc   subroutine corrects for fractional offsets in sampling times
ccc   NOTE : this is not the same as the version of phs_shft used for
ccc   the SEED file verison of dnff

      complex tc,t
      real x(nch,2,nfreq),dr,t1,pi2,dt
      integer nwin
      parameter(pi2 = 6.28318)                      

      do j = 1,nch
ccc      now do the phase shift 
         t1 =   -pi2*dt/(dr*nwin)
         do i = 1,nfreq
            tc = cmplx(cos(i*t1),sin(i*t1))
            t = tc*cmplx(x(j,1,i),x(j,2,i))
            x(j,1,i) = real(t)
            x(j,2,i) = aimag(t)
         enddo
      enddo
      return
      end
