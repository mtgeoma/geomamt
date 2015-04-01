c************************************************************
c
c   this is the subroutine objfun required by the nag routine
c
c     the function is a nonlinear difference
c     equation containing the sums of the differences
c     between the estimated pauli-spin paramters and the
c     measured pauli-spin parameters
c
c     function values are calculated for each frequency
c     (f(*)), the jacobian is also calculated (fjac(*,*)).
c
c     mode    mode of use:
c             mode = 0 only F is computed
c             mode = 1 only FJAC is computed
c             mode = 2 F and FJAC are computed
c
c     neq     no. of equations to solve: 8*(# freqs)*(# sites)
c
c     n       no. of unknowns: 4*(# freqs)*(# sites) + 2*(# sites) + 1
c                                  impedances            ts & sh     azim
c
c     ldfj    first dimension of array FJAC
c
c     x       vector of variables
c
c     f       vector of objective function
c
c     fjac    vector of Jacobian of objective function
c
c     nstate  state
c             nstate = 1 is first call: sets Jacobian to zero
c
c     iuser   workspace integer array
c
c     user    workspace real array
c
c***********************************************************
      subroutine objfun(mode,neq,n,ldfj,x,f,fjac,nstate,
     &      iuser,user)

      implicit none

      include 'size.inc'

      integer n, j, m, p, pe, count(MAXS), nf,
     &        mode, iuser(1), neq, ldfj, i, nstate,
     &        nsite, o

      real*8 th, t, e, c2, s2, ra, user(1),
     &       ia, rb, ib, al_dev(4,MAXF,MAXS),
     &       x((MAXF*4+2)*MAXS+1), f(MAXF*8*MAXS),
     &       fjac(MAXF*8*MAXS,(MAXF*4+2)*MAXS+1), 
     &       dev1, dev2, dev3, dev4

      complex*16 alpha(4,MAXF,MAXS), a1, a2, a3, a4,
     &           za, zb, af1, af2, af3, af4
      complex*16 sigma, delta

      character normtype*6

      logical wantp/.TRUE./, wantq/.TRUE./
      logical DEBUG

      common /cmn1/ alpha, al_dev, normtype
      common /cmn2/ nsite, count
      common /cmndbg/ DEBUG


c---x array
c   1: Re(a) freq. 1 site 1
c   2: Im(a) freq. 1 site 1
c   3: Re(b) freq. 1 site 1
c   4: Im(b) freq. 1 site 1
c   5: Re(a) freq. 2 site 1
c   6: Im(a) freq. 2 site 1
c   7: Re(b) freq. 2 site 1
c   8: Im(b) freq. 2 site 1
c   .
c   .
c   .
c   (no. freqs. * 4) + 2: t for site 1
c   (no. freqs. * 4) + 3: e for site 1
c   (same for other sites)
c   .
c   .
c   .
c   strike for all sites 
    
      if( DEBUG ) then
        write(*,*)'objfun entered: normtype = ', normtype
      endif

c --Set up values for individual difference equations
c   and jacobian

      th = x(n)
      c2 = dcos(2.d0*th)
      s2 = dsin(2.d0*th)

c --loop over sites
  
      p  = 0
      nf = 0
      do o = 1, nsite

        if(o .gt. 1) then
          p = p + count(o-1)*4 + 2
          nf = nf + count(o-1)
        endif
        pe = p + count(o)*4
        t = x(p + count(o)*4 + 1)
        e = x(p + count(o)*4 + 2)

        do m = 1, count(o)

          ra = x(p + (m-1)*4 + 1)
          ia = x(p + (m-1)*4 + 2)
          rb = x(p + (m-1)*4 + 3)
          ib = x(p + (m-1)*4 + 4)
          za = dcmplx(ra,ia)
          zb = dcmplx(rb,ib)
          a1 = alpha(1,m,o)
          a2 = alpha(2,m,o)
          a3 = alpha(3,m,o)
          a4 = alpha(4,m,o)

c...obtain standard deviations (Zxx, Zxy, Zyx, Zyy)
          dev1 = al_dev(1,m,o)
          dev2 = al_dev(2,m,o)
          dev3 = al_dev(3,m,o)
          dev4 = al_dev(4,m,o)

c...define the normalizing factors
          if( normtype(1:2).eq.'L2' ) then
            dev1 = 1.
            dev2 = 1.
          elseif( normtype(1:5).eq.'MAXSD' ) then
            dev1 = dmax1( dev1, dev2, dev3, dev4 )
            dev2 = dev1
          elseif( normtype(1:5).eq.'GAVSD' ) then
            dev1 = ( dev1*dev2*dev3*dev4 )**0.25
            dev2 = dev1
          elseif( normtype(1:6).eq.'GAVSD2' ) then
            dev1 = dsqrt( dev1*dev4 )
            dev2 = dsqrt( dev2*dev3 )
          elseif( normtype(1:5).eq.'SUMSQ' ) then
            dev1 = dsqrt(dev1*dev1 + dev2*dev2 + dev3*dev3 + dev4*dev4)
            dev2 = dev1
          elseif( normtype(1:6).eq.'SUMSQ2' ) then
            dev1 = dsqrt( dev1*dev1 + dev4*dev4 )
            dev2 = dsqrt( dev2*dev2 + dev3*dev3 )
          else
            write(*,*)'Normalizing type not recognised: normtype =',
     &                normtype
            stop
          endif

c...Eqn 35 of GB
	  sigma = za + zb
	  delta = za - zb

c --calculate function values

          if((mode .eq. 0).or.(mode .eq. 2)) then
c...calculate individual differences normalized by the averaged standard 
c   deviations
c

c...expressions for the alphas taken from GB89. 
c...alpha0: eqn 34a (GB)
c           af1 = t*sigma + e*delta
            af1 = - zb*(e-t) + za*(t+e) 

c...alpha1: eqn 34b (GB)
c            af2 = (delta - e*t*sigma)*c2 - (t*delta + e*sigma)*s2
             af2 = s2*(-zb*(e-t))       +c2*(-zb*(1.0d0+t*e))
     &           + c2*( za*(1.0d0-e*t)) -s2*( za*(t+e))
        
c...alpha2: eqn 34c (GB)
c           af3 = -sigma + e*t*delta
            af3 = - zb*(1.0d0+t*e) - za*(1.0d0-e*t)

c...alpha3: eqn 34d (GB)
c           af4 = -(t*delta + e*sigma)*c2 - (delta - e*t*sigma)*s2
            af4 = c2*(-zb*(e-t))       -s2*(-zb*(1.0d0+t*e))
     &           -s2*( za*(1.0d0-e*t)) -c2*(za*(e+t))

            f((nf+m-1)*8 +1) = (dreal(a1)-dreal(af1))/dev1
            f((nf+m-1)*8 +2) = (dimag(a1)-dimag(af1))/dev1
            f((nf+m-1)*8 +3) = (dreal(a2)-dreal(af2))/dev2
            f((nf+m-1)*8 +4) = (dimag(a2)-dimag(af2))/dev2
            f((nf+m-1)*8 +5) = (dreal(a3)-dreal(af3))/dev2
            f((nf+m-1)*8 +6) = (dimag(a3)-dimag(af3))/dev2
            f((nf+m-1)*8 +7) = (dreal(a4)-dreal(af4))/dev1
            f((nf+m-1)*8 +8) = (dimag(a4)-dimag(af4))/dev1

          endif

c...calculate jacobian

          if( (mode.eq.1).or.(mode.eq.2)) then

c...set initial jacobian to zero for first entry

          if( nstate .eq. 1 ) then 
            do i = 1, neq
              do j = 1, n
                fjac(i,j) = 0.0D+00
              enddo
            enddo
            nstate = 2
          endif

c --real a
          fjac((nf+m-1)*8+1,p+(m-1)*4+1) = - (e + t)/dev1
          fjac((nf+m-1)*8+3,p+(m-1)*4+1)=(-c2*(1.0-t*e) + s2*(t+e))/dev2
          fjac((nf+m-1)*8+5,p+(m-1)*4+1) =  (1.0 - e*t)/dev2
          fjac((nf+m-1)*8+7,p+(m-1)*4+1)=(s2*(1.0-e*t) + c2*(t+e))/dev1
c --imaginary a
          fjac((nf+m-1)*8+2,p+(m-1)*4+2) = - (e + t)/dev1
          fjac((nf+m-1)*8+4,p+(m-1)*4+2)=(-c2*(1.0-t*e) + s2*(t+e))/dev2
          fjac((nf+m-1)*8+6,p+(m-1)*4+2) =  (1.0 - e*t)/dev2
          fjac((nf+m-1)*8+8,p+(m-1)*4+2)=(s2*(1.0-e*t) + c2*(t+e))/dev1
c --real b
          fjac((nf+m-1)*8+1,p+(m-1)*4+3) =  (e - t)/dev1            
          fjac((nf+m-1)*8+3,p+(m-1)*4+3)=(s2*(e-t) + c2*(1.0+t*e))/dev2
          fjac((nf+m-1)*8+5,p+(m-1)*4+3) =  (1.0 + e*t)/dev2
          fjac((nf+m-1)*8+7,p+(m-1)*4+3)=(c2*(e-t) - s2*(1.0+e*t))/dev1
c --imaginary b
          fjac((nf+m-1)*8+2,p+(m-1)*4+4) =  (e - t)/dev1            
          fjac((nf+m-1)*8+4,p+(m-1)*4+4)=(s2*(e-t) + c2*(1.0+t*e))/dev2
          fjac((nf+m-1)*8+6,p+(m-1)*4+4) =  (1.0 + e*t)/dev2
          fjac((nf+m-1)*8+8,p+(m-1)*4+4)=(c2*(e-t) - s2*(1.0+e*t))/dev1
c --theta (strike)
          fjac((nf+m-1)*8+3,n) = 2*(c2*(ra*(e+t)+rb*(e-t))
     &                          +s2*(ra*(1.0-t*e)-rb*(1.0+t*e)))/dev2
          fjac((nf+m-1)*8+4,n) = 2*(c2*(ia*(e+t)+ib*(e-t))
     &                          +s2*(ia*(1.0-t*e)-ib*(1.0+t*e)))/dev2
          fjac((nf+m-1)*8+7,n) = 2*(c2*(ra*(1.0-t*e)-rb*(1.0+t*e))
     &                          -s2*(ra*(e+t)+rb*(e-t)))/dev1
          fjac((nf+m-1)*8+8,n) = 2*(c2*(ia*(1.0-t*e)-ib*(1.0+t*e))
     &                          -s2*(ia*(e+t)+ib*(e-t)))/dev1
c --t (twist) 
          fjac((nf+m-1)*8+1,pe+1) = (-rb-ra)/dev1 
          fjac((nf+m-1)*8+2,pe+1) = (-ib-ia)/dev1
          fjac((nf+m-1)*8+3,pe+1) = ((ra-rb)*s2 + e*(ra+rb)*c2)/dev2
          fjac((nf+m-1)*8+4,pe+1) = ((ia-ib)*s2 + e*(ia+ib)*c2)/dev2
          fjac((nf+m-1)*8+5,pe+1) = e*(rb-ra)/dev2
          fjac((nf+m-1)*8+6,pe+1) = e*(ib-ia)/dev2
          fjac((nf+m-1)*8+7,pe+1) = ((ra-rb)*c2 - e*(ra+rb)*s2)/dev1
          fjac((nf+m-1)*8+8,pe+1) = ((ia-ib)*c2 - e*(ia+ib)*s2)/dev1
c --e (shear) 
          fjac((nf+m-1)*8+1,pe+2) = (rb-ra)/dev1 
          fjac((nf+m-1)*8+2,pe+2) = (ib-ia)/dev1
          fjac((nf+m-1)*8+3,pe+2) = ((ra+rb)*s2 + t*(ra+rb)*c2)/dev2
          fjac((nf+m-1)*8+4,pe+2) = ((ia+ib)*s2 + t*(ia+ib)*c2)/dev2
          fjac((nf+m-1)*8+5,pe+2) = t*(rb-ra)/dev2
          fjac((nf+m-1)*8+6,pe+2) = t*(ib-ia)/dev2
          fjac((nf+m-1)*8+7,pe+2) = ((ra+rb)*c2 - t*(ra+rb)*s2)/dev1
          fjac((nf+m-1)*8+8,pe+2) = ((ia+ib)*c2 - t*(ia+ib)*s2)/dev1

        endif

        enddo
      enddo

      return
      end
