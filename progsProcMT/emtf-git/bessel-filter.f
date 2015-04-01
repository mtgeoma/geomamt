      program bessel_filter
      implicit none
      integer n, i, j, nwin
      real t, sf, t0, fac, denom, pi, w
      real q1, q2, q3, w1, w2, w3, ww1, ww2, ww3
      real q(3), ww(3), wwn
      complex cw, resp, resp1, afcor, z(4)

      data q /1.023310,0.611195,0.510318/
      data ww /102.9846,91.32704,86.72125/

      pi=4.*atan(1.)
      sf=0.2
      nwin=65536
      n=nwin/2
      t0=0.2

      ww1 = 102.9846
      ww2 = 91.32704
      ww3 = 86.72125
      q1 = 1.023310
      q2 = 0.611195
      q3 = 0.510318

      fac = 15.0 * t0
      
      w1 = ww1/fac
      w2 = ww2/fac
      w3 = ww3/fac
      z(1)=(w1*w1)*(w2*w2)*(w3*w3)
      do i=1,n
         t=nwin/(i*sf)
         w=2*pi/t
         cw=cmplx(0.0,w)
         z(2)=1.0/(cw**2 + cw*w1/q1 + w1**2)
         z(3)=1.0/(cw**2 + cw*w2/q2 + w2**2)
         z(4)=1.0/(cw**2 + cw*w3/q3 + w3**2)
         resp = z(1)*z(2)*z(3)*z(4)
         denom = CABS(resp)**2
         resp=resp/denom

         afcor=cmplx(1.0,0.0)
         do j=1,3
            wwn=ww(j)/fac
            afcor = afcor*wwn**2/(cmplx(wwn**2-w**2,w*wwn/q(j)))
         enddo
         denom=CABS(afcor)**2
         afcor = denom/afcor
C         write(*,*) afcor,1.0/resp
         write(*,*) real(afcor),imag(afcor),real(1.0/resp),
     +              imag(1.0/resp)
      enddo
      end
