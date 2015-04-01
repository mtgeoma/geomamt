c*************************************
        function afcor(iftype,t0,t,alpha)
        
c       computes filter responses for four types of filters used in
c       emslab long period mt experiment

cf      7.6.95  LMT MAGSON MAGNETOMETER low pass (10 sec) added

c       iftype  = 1 ==> no filter
c               = 2 ==> 1-pole high pass (uw mt)
c               = 3 ==> 2-pole high pass (epb mt)
c               = 4 ==> 2-pole low pass (uw,epb anti alias)
c               = 5 ==> 1-pole low pass (mag integrators)
c               = 6 ==> "box car" low pass filter - used by filloux
C-PDAS          =11 ==> PDAS 200 Hz low pass
C-PDAS          =12 ==> GMS05 coil system function
cf              =13 ==> LMT RAP Bessel low pass (10 sec) 
cf              =14 ==> LMT MAGSON MAGNETOMETER low pass (5 sec)
cf              =15 ==> LMT RAP Bessel low pass (120 sec)
cf              =16 ==> LMT MAGSON MAGNETOMETER low pass (10 sec)
c               =19 ==> LRMT Bessel filter

        complex afcor,temp
C-PDAS  LMT
        complex ce, p, pl
        real butter(3), butt2(4)
        real creal, cimag, acoef(10), bcoef(10), lmta(4), lmtb(5)
        real m1acoef(10),m1bcoef(10)
        real lmt120a(4), lmt120b(5)

c       LRMT Bessel filter
        real q(3), ww(3), wwn, w, pi, denom
        data q /1.023310,0.611195,0.510318/
        data ww /102.9846,91.32704,86.72125/

        data lmta/0.9995,0.4884,-540.44,3314.08/
        data lmtb/0.0003,-30.2821,131.976,-645.841,56054.2/

        data butt2/1.,1.8019,1.2470,0.4450/
        data butter /1.9318,1.4141,0.5176/

cf......LMT RAP Bessel low pass (120 sec)
        data lmt120a/1.00938, -17.4907,-72672.1, 5.9892e06/
        data lmt120b/-0.0003034, -363.385, 19003.1, -1.115e06, 1.162e09/

cf......LMT MAGSON MAGNETOMETER low pass 
C        *********************************
c        5sec:
c        *********************************

      data acoef/0.9999288296,0.0698487366,-36.54067611,
     $       -48.07610115,1362.242668,-5500.831271,10839.27048,
     $       -11784.68695,6791.730771,-1624.96809/
      data bcoef/-0.001502197016,5.873206596,46.9447023,
     $       -705.7354187,3305.132341,-8326.369934,12573.32541,
     $       -11432.56764,5792.333082,-1259.449763/

c        data acoef/5.1246, 0.86734, -522.35, 2967.31, -5167.64/
c        data bcoef/0.03539, -52.1911, 81.38222, 1594.21, -5209.82/

c        *********************************
c         10sec:
c        *********************************
      data m1acoef/0.9997296956,0.6799243865,-194.6183338,
     $       1829.13399,-8241.197348,21653.75418,-34803.94742,
     $       33677.68098,-17990.74652,4066.264465/

      data m1bcoef/-0.02305288039,17.50142214,-141.2599948,
     $       268.8431802,1396.30022,-9193.972813,23246.24695,
     $       -31054.26697,21651.84395,-6217.244816/

c        data m1acoef/0.999608, 0.784654, -181.402, 1602.82,
c     &               -6747.79, 16490.1, -24583.9, 22046.0,
c     &               -10927.5 ,2299.88/
c        data m1bcoef/-0.0167768,15.7718,-112.219,148.612,1104.43,
c     &               -5478.73,11206.8,-12160.5,6864.87,-1589.04/

        ce = cmplx(1.,0.)

C-PDAS END
        go to(100,200,300,400,500,600,700,800,900,1000,1100,1200,
     &        1300,1400,1500,1600,1700,1800,1900) iftype

C-PDAS
700     return
800     return
900     return
1000    return
1700    return
1800    return
C-PDAS END
        
100     continue
c       no filter
        afcor = (1.0,0.0) 
        return

200     continue
c               uw-mt high pass
        gain = 1./sqrt(1.+(t/t0)*(t/t0))
        phase = atan(t/t0)
        temp = cmplx(0.,phase)
        afcor = gain*cexp(temp)
        return
        
300     continue
        
        r=t/t0
        t1 = r**4. + (alpha**2.-2.)*r**2. + 1.
        afcor = cmplx( (1.-r*r)/t1, alpha*r/t1 )
        return

         
400     continue
        
        r=t0/t
        t1 = r**4. + (alpha**2.-2.)*r**2. + 1.
        afcor = cmplx( (1.-r*r)/t1, (-1.)*alpha*r/t1 )
        return
        
500     continue
        r = t0/t
        gain = 1./sqrt(1.+(r*r))
        phase = atan(-r)
        temp = cmplx(0.,phase)
        afcor = gain*cexp(temp)
        return

600     continue
c    here t0 is width of box-car, alpha gives offset of boxcar center
c       from nominal sampling time
        if (t0.eq.0.) then
           afcor = 1.0
        else
           r = (t0/t)*3.14159
           afcor = sin(r)/r
        end if
                             
        temp = cmplx(0.,6.28318*alpha/t)
        afcor = afcor * cexp(temp)
        return

C-PDAS  PDAS-Butterworth filter

1100    continue
        afcor = (1.0,0.0) 
        r= 1./(t*200)
        do 1101 ibw = 1,3
1101    afcor = afcor*(cmplx(1.-r*r,butter(ibw)*r))
        return

C-PDAS  GMS05 coil sys. function
1200    p = cmplx(0.,1./(t*4.))
        pl = cmplx (0.,1./(t*8192.))
        afcor = (800.*(p/(ce+p))*(ce/(ce+pl)))
c       print*,t, 1./t, cabs(afcor), atan2(aimag(afcor),real(afcor))*180./3.1415
c       ce/
        return
C-PDAS END


        
C LMT



C       RA 7-pol Butterworth lowpass at 10 sec
C1300   r = 10./t
C       afcor = cmplx ( 1.  ,     butt2(1) *  r  )
C    &        * cmplx ( 1. - r*r, butt2(2) *  r  )
C     &        * cmplx ( 1. - r*r, butt2(3) *  r  ) 
C     &        * cmplx ( 1. - r*r, butt2(4) *  r  )
C       afcor = 1./afcor
c       print*,'T, abs, phi',t,(cabs(afcor)),(180./3.14159*
c     &  atan(aimag(afcor)/real(afcor)))



C       return


C       Fitting RA 7-Pol Bessel filter by a 4 and 3 pol polynom
1300    creal = 0.
        cimag = 0.
        do 1301 i = 1,4
        creal = creal + lmta(i)*(1./t)**(i-1)        
1301    continue
        do 1302 i = 1,5
        cimag = cimag + lmtb(i)*(1./t)**(i-1)
1302    continue
        afcor =  (cmplx(creal,cimag))
c       print*,'T, abs, phi',t,(cabs(afcor)),(180./3.14159*
c     &  atan(aimag(afcor)/real(afcor)))


        return


C      MA  MAGSON magnetometer transfer function

1400    creal = 0.
        cimag = 0.
        do 1411 i = 1,10
        creal = creal + acoef(i)*(1./t)**(i-1)
        cimag = cimag + bcoef(i)*(1./t)**(i-1)
1411    continue
        afcor = (cmplx(-creal,cimag))/25.
c       print*,'T, abs, phi',t,(cabs(afcor)),(180./3.14159*
c     &  atan(aimag(afcor)/real(afcor)))


        return


C       R2 RAP 7 - Pol Bessel Low Pass at 120 s

1500    creal = 0.
        cimag = 0.
        do 1501 i = 1,4
        creal = creal + lmt120a(i)*(1./t)**(i-1)        
1501    continue
        do 1502 i = 1,5
        cimag = cimag + lmt120b(i)*(1./t)**(i-1)
1502    continue
        afcor =  (cmplx(creal,cimag))
c       print*,'RA: T, abs, phi',t,(cabs(afcor)),(180./3.14159*
c     &  atan(aimag(afcor)/real(afcor)))

        return

C       M1 MAGSON magnetometer transfer function (10 sec low pass)

1600    creal = 0.
        cimag = 0.
        do 1601 i = 1,10
        creal = creal + m1acoef(i)*(1./t)**(i-1)
        cimag = cimag + m1bcoef(i)*(1./t)**(i-1)
1601    continue
        afcor = (cmplx(creal,cimag))
c       print*,'M1: T, abs, phi',t,(cabs(afcor)),(180./3.14159*
c    &  atan(aimag(afcor)/real(afcor)))

        return

 1900   pi=4.*atan(1.)
        w=2*pi/t
        afcor=cmplx(1.0,0.0)
        do i = 1,3
           wwn=ww(i)/(15.0*t0)
           afcor = afcor*wwn**2/(cmplx(wwn**2-w**2,w*wwn/q(i)))
        enddo
        denom=CABS(afcor)**2
        afcor = denom/afcor

        return

        end
