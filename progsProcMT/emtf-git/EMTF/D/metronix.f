      subroutine startmtx(mode,nfi,ffreq,freal,fimag,yp1,yp2,p)
      include 'iounits.inc'
      include 'metronix.inc'
      integer nfi,ierr,i,i1,i2,j1,j2,k1,k2
      character par1*80, par2*80,mode*3
      real ffreq(nfilm),ffreq1(nfilm),freal(nfilm),fimag(nfilm),
     $     yp1(nfilm),yp2(nfilm),temp(nfilm)
      real amp,phi,pi,rad,p(3)
      complex p1c,p2c,respmtx
      parameter (pi=3.141592653589793,rad=pi/180.)

      call strtolower(mode)
      call strlen(cfrsp,i1,i2)
      open(unit=rsp_unit,file=cfrsp(i1:i2),status='old')

      call strlen(mode,i1,i2)
      ! do until find a line that start with 'chopper' and 'mode'
 4    continue
      read(rsp_unit,*)par1,par2
      call strtolower(par1)
      call strtolower(par2)
      call strlen(par1,j1,j2)
      call strlen(par2,k1,k2)
      if (par1(j1:j2).ne.'chopper'.or.par2(k1:k2).ne.mode(i1:i2)) goto 4

      ! read filter frequency, amplitude (normalized by freq) and phase
      ! until end of file, an error occurs (didn't find 3 floats)
      ! or reach maximum number of filter parameters (nfilm)
      do i=1,nfilm
         read(rsp_unit,*,end=6,err=6)ffreq(i),amp,phi
         freal(i)=amp*cos(rad*phi)*ffreq(i)
         fimag(i)=amp*sin(rad*phi)*ffreq(i)
      enddo
 6    nfi=i-1
      close(unit=rsp_unit)

      ! include two extra points from theoretical transfer function
      ! to stabilize interpolation at high frequency
      amp=1.5*ffreq(nfi)/ffreq(nfi-1)
      do i=1,2
         if(nfi.eq.nfilm) then
            exit
         else
            nfi=nfi+1
            ffreq(nfi)=ffreq(nfi-1)*amp
            p1c = cmplx(0.,ffreq(nfi)/p(2))
            p2c = cmplx(0.,ffreq(nfi)/p(3))
            respmtx = p(1)*(p1c/(1.0+p1c))*(1.0/(1.0+p2c))
            freal(nfi)=real(respmtx)
            fimag(nfi)=imag(respmtx)
         endif
      enddo

      do i=1,nfi
         ffreq1(i)=ffreq(i)
      enddo
      call ssort(ffreq,freal,nfi,2)
      call ssort(ffreq1,fimag,nfi,2)
      call curv1(nfi,ffreq,freal,0,0,3,yp1,temp,sigma,ierr)
      if(ierr.ne.0)then
         write(*,*)'error in first curv1:',ierr
         stop
      endif
      call curv1(nfi,ffreq,fimag,0,0,3,yp2,temp,sigma,ierr)
      if(ierr.ne.0)then
         write(*,*)'error in second curv1:',ierr
         stop
      endif
      end

      function respmtx(nfi,ffreq,freal,fimag,yp1,yp2,p,
     $     mode,fr)
      include 'iounits.inc'
      include 'metronix.inc'
      integer nfi,i1,i2
      real r,q,p(3)
      character mode*3
      complex p1c,p2c,respmtx
      real ffreq(nfilm),freal(nfilm),fimag(nfilm),yp1(nfilm),
     $     yp2(nfilm)

      call strtolower(mode)
      call strlen(mode,i1,i2)
      if( (mode(i1:i2).eq.'ttf') .or.
     $    (mode(i1:i2).eq.'on'.and.fr<ffreq(1)) ) then
         ! use theoretical transfer function
         p1c = cmplx(0.,fr/p(2))
         p2c = cmplx(0.,fr/p(3))
         respmtx = p(1)*(p1c/(1.0+p1c))*(1.0/(1.0+p2c))*1000.0
      else
         r=curv2(fr,nfi,ffreq,freal,yp1,sigma)*1000.0
         q=curv2(fr,nfi,ffreq,fimag,yp2,sigma)*1000.0
         respmtx=cmplx(r,q)
      endif
      return
      end
