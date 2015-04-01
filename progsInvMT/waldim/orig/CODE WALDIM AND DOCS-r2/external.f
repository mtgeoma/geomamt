C   *****************************************************************

C    External functions and subroutines

C   *****************************************************************

C     Computing invariants I1, I2, I3, I4

      double precision function Yinv12(a,b)
      implicit double precision (a-h,o-z)
        Yinv12=(a**2+b**2)**0.5   
      return
      end

      double precision function Yinv34(a,b,c)
      implicit double precision (a-h,o-z)
        Yinv34=((a**2+b**2)**0.5)/c  
      return
      end

C   *****************************************************************
C     First quadrant conversion

      double precision function FQ(a)
      implicit double precision (a-h,o-z)
      
   10 if(a.gt.90) then
      a=a-90
      else
      goto 20
      endif      
      goto 10
      
   20 if(a.lt.0) then
      a=a+90
      else
      goto 30
      endif
      goto 20
   30 fq=a
      return 
      end
           
C   *****************************************************************

c Program to produce exponentially correlated colored (Gaussian) noise.
c based on Fox et al Physical Review A vol.38(1988)5938 and 
      
      double precision function gasdev(idum)
c   function gasdev
c   (Numerical Recipes, W.T.Vetterling, et.al., p.203)
c
c   description
c   ===========
c   returns a normally distributed deviate with zero mean and unit variance
c   using usran(idum) in file ranims.f
c   (formerly used ran1(idum) in file ran1ol.f)
c
c 09/01/88:converted to double precision
c 07/06/89:added save for gset (most computers dont need these save statements)
      implicit double precision (a-h,o-z)

      save iset, gset
      data iset/0/
      if (iset.eq.0) then
c***we don't have an extra deviate handy,
c***so pick two uniform random numbers in the square extending from -1 to +1
c***in each direction.
1       v1=2.d0*usran(idum)-1.d0
        v2=2.d0*usran(idum)-1.d0
        r=v1**2+v2**2
        if(r.ge.1.d0)go to 1
        fac=sqrt(-2.d0*log(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
c***we do have an extra deviate handy,
c***so return it, and unset the flag.
        gasdev=gset
        iset=0
      endif
      return
      end


C   *****************************************************************

      double precision function usran(ir)
c
c   this subroutine generates random values between 0.0 and 1.0 using
c   an integer seed
c   it is based on the imsl routine ggubs.
c
c   double precision version
c
      implicit double precision (a-h,o-z)
      parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
      iir=abs(mod(da*ir,db)+0.5d0)
      ir = iir
      usran=float(ir)/dc
      
      return
      end


C   *****************************************************************

      subroutine nullfreq (mp,nfreq,f,rot,zr,zi,zvar)

      implicit double precision (a-h,o-z)
      Dimension f(mp),rot(mp),zr(2,2,mp),zi(2,2,mp),zvar(2,2,mp)
      integer nfreq,nfreq2
 
      pnull1=1.00e31

      nfreq1=nfreq

      do k=1,nfreq1

      suma=0
      do i=1,2
        do j=1,2
        suma=suma+abs(zr(i,j,k))+abs(zi(i,j,k))+abs(zvar(i,j,k))
        end do
      end do
      
      suma=suma + abs(f(k))+abs(rot(k))
           
        if (k.le.nfreq) then

        if (suma.ge.pnull1) then
           nfreq2=nfreq-1

           do k2=k,nfreq2
                   
         f(k2)=f(k2+1) 
       rot(k2)=rot(k2+1)
       do i=1,2
         do j=1,2
                 zr(i,j,k2)=zr(i,j,k2+1)
                 zi(i,j,k2)=zi(i,j,k2+1)
                zvar(i,j,k2)=zvar(i,j,k2+1)
           end do
             end do


           end do
          f(nfreq)=0 
          rot(nfreq)=0
            do i=1,2
              do j=1,2
                 zr(i,j,nfreq)=0
                 zi(i,j,nfreq)=0
                 zvar(i,j,nfreq)=0
            end do
             end do

        nfreq=nfreq2
        end if
        else
        goto 40
        end if

        end do

 40     return
      end
