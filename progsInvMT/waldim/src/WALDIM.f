C***************************************************************************************
C     WALDIM
C
C     A. Marti, P. Queralt and J. Ledo, 2008, "WALDIM: A code for the dimensionality 
C     analysis of magnetotelluric data using the Rotational Invariants of the 
C     Magnetotelluric Tensor",
C     Submitted to Computers and Geosciences.
C
C   ***************************************************************************************
C
C     Input Files:
C     
C     1:param.cfg: configuration file, if it does not exists, the program
C     creates it.
C
C     11: list with sites and coordinates
C
C     2: Data Files *00i.EDI with impedances (one EDI file for each site)
C  
C     Output files:
C
C     3:  inv.dat: table with the invariant true values
C     4:  cases***.dat: file with dimensionality for every site and period
C
C     7: err***.dat: invariants, strikes and errors
C     8: file with true values, stat. values and errors, and biases (optional)
C     9: dimensionality results for each period and band classification for each site

C     10: file with other dimensionality indicators

C     14: dimensionality classification, all bands
C     15: dimensionality classification, one file per band
C

C
C   ***************************************************************************************
C
C   Relevant parameters and data description:
C   -----------------------------------------
C
C   *** Character strings: ***
C
C       realt(3): character vector containing date, time and zone.
C
C       File11,File2,File3,File4,File7,File8,File9,File10,File14,File15:Char.strings to identify file names 
C
C     sidsite: Char.strings to identify site names
C
C     eif, uef, avg, Units, statf, output: input options
C
C     str: auxiliary strings 
C
C     Info: long character strings to be read in files
C
C     namedim(8)*9: char. string to identify dimensionality types
C
C
C   *** Integers: ***
C
C     icas: dimensionality cases, idenfified from 0 to 7:
C         0: undetermined, 
C         1: 1D
C         2: 2D
C         3: 3D/2D twist
C         4: 3D/2D
C         5: 3D 
C         6: 3D/2Ddiag 
C         7: 3D/1D2D
C
C nsites,nsit: nr of sites
C
C   *** Parameters: ***
C
C     mp: max index in period dimension arrays (default mp=100)
C     mp: max index in site dimension arrays (default ms=100)
C     mnrel: max number of realizations in statistical realizations (default mnrel = 1000)
C
C   *** Real: ***
C     implicit double precision (a-h,o-z)
C
C     Cts: constant used for unit conversions
C
C     df, dm: auxiliary functions to compute derivatives
C     th, thq: thresholds for invariants I3-I7 and Q.
C     pm, pma, bpd: periods min., max. and nr of bands for every decade in the average
C
C   *** Vectors and matrices ***
C
C     jdate: vector containing data and time 
C
C     nf: number of frequencies at each site
C
C     lat,long, plati, plongi: latitudes and longitudes
C
C     f,af, per: frequency and period values at each site
C
C     rot,arot, zr, zi, zvar: rotation, real and imaginary impedances and variances of the data
C
C     errel: mean relative error of the impedance components
C
C     fi, pnu, di: linear combinations of the impedance components
C
C     efn: errors of fi and pnu
C
C     derf, derm: derivatives of fi and pnu used to compute errors
C
C     Y:  invariants
C
C     Yst, aYst: parameters related to invariants
C
C     ry, rstr: statistical realizations of invariants and related parameters
C
C     DevY, devst,adevst: errors of invariants and related parameters
C
C     Bias, biasst: biases between true and statistical values of invariants and related params.
C
C     yc: invariant classification as 0, 1, undetermined (1000) or unclassified (34000)
C
C     nicas: number of cases of each dimensionality type
C
C     aux: auxiliary vectors to be read
C
C     iacas: dimensionality type of each frequency
C
C     ict: dimensionality mode in a frequency band
C
C     *** External programs: ***
C
C     external.f: contains external functions and subroutines
C     inoutdata.f: contains subroutines used for input/output 
C
C     *** Functions: *** (external.f)
C
C     Yinv12, Yinv34: compute invariants 1,2, 3 and 4.
C     FQ: converts angles into the first quadrant
C     Gasdev: generates random gaussian noise
C     Usran: function used in function Gasdev
C     NullFreq: searches for null values in edi files (default: NULL= 1E32)
C 
C     *** Subroutines: ***  (inoutdata.f)
C     
C     readparams: reads parameters from param.cfg file or from keyboard

C     readtbl: reads table with sites and coordinates 
C     Outputheaders: generates headers of output files
C     Header14: generates header of file14
C     Header15: generates header of file15
C  
C     writefile3,4,7,8,14: writes results in output files
C
C     Header9 and Header9bis: generates headers of file9
C     WriteFile9 and WriteFile9bis: writes results in file9
C
C     CompWriteFile10: computes and writes parameters in file 10
C
C
C
C     *** Intrinsic libraries: ***

C     DATE_AND_TIME (FORTRAN 95/90 INTRINSIC LIBRARY)
C
C   *******************************************************************

C Declaration of variables

      implicit double precision(a-h,o-z)

      Parameter (mp=100,ms=100, mnrel=1000)

      Integer icas,iacas(mp,ms),ict(7)

      Character output*25,sidsite(ms)*25,sidsitet*25
      Integer nout,lsidsite(ms),lsidsitet

      Character realt(3)*12

      Character file11*30,file2*30
      Character file9(ms)*30
      Character namedim(8)*9

      Character str*3
      Character eif*1,uef*1,avg*1,statf*1,Units*1
      Character Info*100

      Dimension jdate(8)
      Dimension nf(ms),plati(ms),plongi(ms)
      Dimension F(mp),Per(mp),Rot(mp),Zr(2,2,mp),Zi(2,2,mp),Zvar(2,2,mp)
      Dimension errel(mp),efn(4,mp)
      Dimension Fi(4,mp),Pnu(4,mp),Di(4,4,mp)
      Dimension derf(4,4,4),derm(4,4,4)
      Dimension Y(10,mp,2),Yst(10,mp,2),ry(10,mnrel),rstr(10,mnrel)
      Dimension Sumy(10),dify(10),sumst(10),difst(10)
      Dimension DevY(10,mp,2),devst(10,mp,2)
      Dimension Bias(8),biasst(10)

      Dimension af(mp,ms),arot(mp,ms),ayst(9,mp,ms)
      Dimension adevst(9,mp,ms)

      Dimension yc(8,mp),nicas(0:7)



C Parameters Pi and magnetic permeability

      Pi=acos(-1.)
      pMu0=4*Pi*1e-7

C Default number of realizations in statistical computations
      nrel=100


      write(6,*) ' **************************************************'
      write(6,*) ' *                   WALDIM                       *'
      write(6,*) ' *   "A code for the dimensionality analysis      *'
      write(6,*) ' *         of magnetotelluric data using          *'
      write(6,*) ' *       the Rotational Invariants                *'
      write(6,*) ' *      of the Magnetotelluric Tensor"            *'
      write(6,*) ' *    by A. Marti, P. Queralt and J. Ledo         *'
      write(6,*) ' *    Computers and Geosciences (submitted 2008)  *'
      write(6,*) ' *                                                *'
      write(6,*) ' **************************************************'
 

C Reading input files and input parameters.

      CALL ReadParams(file11,units,th,thq,eif,uef,ep,
     +               avg,pm,pma,bpd,output,statf,errcfg)

      if (errcfg.eq.1) goto 700

        do i=1,25
          if (output(i:i).eq.' ') then
             nout=i-1
             goto 19
          end if
        end do

C Date and time: seed for random numbers generation
               
  19  CALL date_and_time(realt(1),realt(2),realt(3),jdate)    


C   ***************************************************************************************

C MAIN LOOP FOR EACH SITE: READ EDI FILES.
C FOR EACH FREQUENCY: COMPUTE INVARIANTS, RELATED PARAMETERS, AND DIMENSIONALITY.

C Reading list of sites. 
C     1: info and number of sites

      open(unit=11,file=file11,status='old')
      read(11,'(a)') info
      read(11,*) nsites

C     2: reading site identifications, and cordinates

      lmiss=0
      l=0
      DO 20 le=1,nsites
        read(11,'(a)',end=101) info(1:100)
        ispace=0
        
        do kc=1,100
          if (info(kc:kc).eq.' ') ispace=ispace+1
        end do

        if (ispace.eq.100) goto 101

        CALL readtbl(info,p1,p2,sidsitet,lsidsitet)!READS IDSITE, LAT, LONG

      
C Reading .edi file for site le:
      
        file2=sidsitet(1:lsidsitet)//'.edi'
     
        open(unit=2,file=file2,status='old',err=201)

        l=l+1
        plati(l) = p1
        plongi(l) = p2
        sidsite(l)=sidsitet(1:lsidsitet)
        lsidsite(l)=lsidsitet
   
C Opening and creating headers of output files (only once)
      
          if (l.eq.1) then
      
             CALL Outputheaders(output,nout,th,statf,
     +                    errov)

           if (errov.eq.1) goto 700
        end if

        do 30 i=1,1500
          read(2,'(a)',err=203) Info
          do j=1,10
            if (info(j:j+5).eq.'NFREQ=') goto 31
          end do
   30   continue
      
   31   backspace(2) 

        do ja=j+6,j+10
          if ((info(ja:ja).ne.' ')) goto 32          
        end do

   32   do je=ja,ja+3
          if ((info(je:je).eq.' ')) goto 33
        end do

   33   str=info(ja:je-1)
      
        read (str,'(i3)') nfreq
         
        do 36 i=1,200
          read(2,'(a)',err=203) Info
          do j=1,10
            if (info(j:j+5).eq.'>FREQ') goto 37
          end do
   36   continue
      
   37   Nf(l)=nfreq
    
        read(2,*) (F(k),k=1,Nfreq)

   38   read(2,'(a)') info

        do j=1,10
          if (info(j:j+4).eq.'>ZROT') then
             read(2,*) (Rot(k),k=1,Nfreq)
             goto 38
           end if
        end do

        do j=1,10
          if ((info(j:j+4).eq.'>ZXXR').or.(info(j:j+4).eq.'>Zxxr')) then
             read(2,*) (Zr(1,1,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+4).eq.'>ZXXI').or.(info(j:j+4).eq.'>Zxxi')) then
             read(2,*) (Zi(1,1,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+7).eq.'>ZXX.VAR'). 
     +    or.(info(j:j+7).eq.'>Zxx.var')) then
             read(2,*) (Zvar(1,1,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+4).eq.'>ZXYR').or.(info(j:j+4).eq.'>Zxyr')) then
             read(2,*) (Zr(1,2,k),k=1,Nfreq)
             goto 38
            end if
        end do

        do j=1,10
          if ((info(j:j+4).eq.'>ZXYI').or.(info(j:j+4).eq.'>Zxyi')) then
             read(2,*) (Zi(1,2,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+7).eq.'>ZXY.VAR').
     +    or.(info(j:j+7).eq.'>Zxy.var')) then
             read(2,*) (Zvar(1,2,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+4).eq.'>ZYXR').or.(info(j:j+4).eq.'>ZYxr')) then
             read(2,*) (Zr(2,1,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+4).eq.'>ZYXI').or.(info(j:j+4).eq.'>Zyxi')) then
             read(2,*) (Zi(2,1,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+7).eq.'>ZYX.VAR'). 
     +    or.(info(j:j+7).eq.'>Zyx.var')) then
             read(2,*) (Zvar(2,1,k),k=1,Nfreq)
             goto 38
          end if
        end do

        do j=1,10
          if ((info(j:j+4).eq.'>ZYYR').or.(info(j:j+4).eq.'>Zyyr')) then
             read(2,*) (Zr(2,2,k),k=1,Nfreq)
             goto 38
            end if
        end do
 
        do j=1,10
          if ((info(j:j+4).eq.'>ZYYI').or.(info(j:j+4).eq.'>Zyyi')) then
             read(2,*) (Zi(2,2,k),k=1,Nfreq)
             if ((eif.eq.'n').or.(eif.eq.'N')) then
                goto 404
             else
                goto 38
             end if
          end if
        end do

        do j=1,10
          if ((info(j:j+7).eq.'>ZYY.VAR').or.
     +       (info(j:j+7).eq.'>Zyy.var')) then
             read(2,*) (Zvar(2,2,k),k=1,Nfreq)
             goto 404
          end if
        end do

        goto 38
 
 404      do 40 i=1,2
          do 45 j=1,2
            do k=1,nfreq
              if (Zr(i,j,k).eq.0) then
                 Zr(i,j,k)=1e-6
              end if
              if(Zi(i,j,k).eq.0) then
                 Zi(i,j,k)=1e-6
              end if
      
                if ((uef.eq.'n').or.(uef.eq.'N')) then
                 Zvar(i,j,k)=(ep*0.01*max(Zr(i,j,k),Zi(i,j,k)))**2
              end if
            end do
  45      continue
  40    continue
        close(2)

C Detecting null frequencies:

      Call NullFreq(mp,nfreq,f,rot,zr,zi,zvar)

      nf(l)=nfreq


  
C     Mean relative error at each frequency:

        do 50 k=1,nfreq
          errel(k)=0
          do 55 i=1,2
            do 56 j=1,2
              errel(k)=errel(k)+max(sqrt(abs(Zvar(i,j,k)/Zr(i,j,k))),
     +        sqrt(abs(Zvar(i,j,k)/Zi(i,j,k))))
   56       continue
   55     continue  
          errel(k)=errel(k)*100/4  
   50   continue  
    

C Writing results for each site
        CALL Header9(sidsite(l),lsidsite(l),th,file9(l))


C     Theoretical computation of the invariants, derivatives and errors,
C     using classical error propagation:

C     Units:
C     Cts=1/1000 if units in the impedances are m/s  (M: magnetotelluric tensor)
C     Cts=1 if units in the impedances are km/s (F: Field units)
C     Cts=1/(1000*Mu0) for impedances, (Z: impedance tensor, Z=Mu0*M)

 
        if ((units.eq.'m').or.(units.eq.'M')) then
           Cts=1/1000.
        end if
        if ((units.eq.'f').or.(units.eq.'F')) then
           Cts=1.
        end if
        if ((units.eq.'z').or.(units.eq.'Z')) then
           Cts=1./(1000*pMu0)
        end if
  
C Loop to compute invariant parameters and errors for each frequency
 
        do 60 k=1,nfreq

          lc=lc+1
          Fi(1,k)=0.5*(Zr(1,1,k)+Zr(2,2,k))*Cts
          Pnu(1,k)=0.5*(Zi(1,1,k)+Zi(2,2,k))*Cts
          Fi(2,k)=0.5*(Zr(1,2,k)+Zr(2,1,k))*Cts
          Pnu(2,k)=0.5*(Zi(1,2,k)+Zi(2,1,k))*Cts
          Fi(3,k)=0.5*(Zr(1,1,k)-Zr(2,2,k))*Cts
          Pnu(3,k)=0.5*(Zi(1,1,k)-Zi(2,2,k))*Cts
          Fi(4,k)=0.5*(Zr(1,2,k)-Zr(2,1,k))*Cts
          Pnu(4,k)=0.5*(Zi(1,2,k)-Zi(2,1,k))*Cts
          
C Errors:
          
          do 61 ld=1,4
            if ((ld.eq.1).or.(ld.eq.3)) then
               efn(ld,k)=0.5*sqrt((Zvar(1,1,k))+(Zvar(2,2,k)))*Cts
            else
               efn(ld,k)=0.5*sqrt((Zvar(1,2,k))+(Zvar(2,1,k)))*Cts
            end if
   61     continue

C Invariants and errors using classical error propagation (index:1):
     
          Y(1,k,1)=Yinv12(fi(4,k),fi(1,k))
          Y(2,k,1)=Yinv12(pnu(4,k),pnu(1,k))
          Y(3,k,1)=Yinv34(fi(2,k),fi(3,k),Y(1,k,1))
          Y(4,k,1)=Yinv34(pnu(2,k),pnu(3,k),Y(2,k,1))

C Parameters dij (differences between fis and pnus: 
          do 62 m1=1,4
            do 63 m2=1,4
              if (m2.ne.m1) then
                di(m1,m2,k)=(fi(m1,k)*pnu(m2,k)-pnu(m1,k)*fi(m2,k))/
     +          (Y(1,k,1)*Y(2,k,1))
             
C These derivaives will be useful for d12,d13,d34 i d24.

                do 64 m3=1,4
                  if (m3.eq.m1) then
                     df1=pnu(m2,k)
                     dm1=-fi(m2,k)
                  end if
                  if (m3.eq.m2) then
                     df1=-pnu(m1,k)
                     dm1=fi(m1,k)
                  end if
                  if ((m3.ne.m1).and.(m3.ne.m2)) then
                     df1=0
                     dm1=0
                  end if
                  if ((m3.eq.1).or.(m3.eq.4)) then
                     df2=fi(m3,k)/Y(1,k,1)
                     dm2=pnu(m3,k)/Y(2,k,1)
                  else
                     df2=0
                     dm2=0
                  end if
                  derf(m1,m2,m3)=df1/(Y(1,k,1)*Y(2,k,1))-
     +            di(m1,m2,k)*df2/Y(1,k,1)
                  derm(m1,m2,m3)=dm1/(Y(1,k,1)*Y(2,k,1))-
     +            di(m1,m2,k)*dm2/Y(2,k,1)
                
   64           continue
              end if
   63       continue
   62     continue

          Y(5,k,1)=(fi(4,k)*pnu(1,k)+fi(1,k)*pnu(4,k))
     +    /(Y(1,k,1)*Y(2,k,1))
          Y(6,k,1)=di(4,1,k)
          Y(8,k,1)=sqrt((di(1,2,k)-di(3,4,k))**2+
     +    (di(1,3,k)+di(2,4,k))**2)
          
          do j=1,8
            if (y(j,k,1).eq.0) then
               y(j,k,1)=1e-5
            end if
          end do

          Y(7,k,1)=(di(4,1,k)-di(2,3,k))/Y(8,k,1)
          
          devY(1,k,1)=sqrt((fi(1,k))**2*(efn(1,k))**2+(fi(4,k))**2*
     +    (efn(4,k))**2)/Y(1,k,1)
          
          devY(2,k,1)=sqrt((pnu(1,k))**2*(efn(1,k))**2+(pnu(4,k))**2*
     +    (efn(4,k))**2)/Y(2,k,1)
          
          devY(3,k,1)=sqrt((1/(Y(3,k,1))**2)*((fi(2,k))**2*
     +    (efn(2,k))**2+(fi(3,k))**2*(efn(3,k))**2)+((Y(3,k,1))**2)*
     +    ((fi(1,k))**2*(efn(1,k))**2+(fi(4,k))**2*(efn(4,k))**2))
     +    /Y(1,k,1)**2
          
          devY(4,k,1)=sqrt((1/(Y(4,k,1))**2)*((pnu(2,k))**2*
     +    (efn(2,k))**2+(pnu(3,k))**2*(efn(3,k))**2)+((Y(4,k,1))**2)*
     +    ((pnu(1,k))**2*(efn(1,k))**2+(pnu(4,k))**2*(efn(4,k))**2))
     +    /Y(2,k,1)**2
          
          devY(5,k,1)=(1/(Y(1,k,1)*Y(2,k,1)))*
     +    sqrt((pnu(4,k)-Y(5,k,1)*Y(2,k,1)*fi(1,k)/Y(1,k,1))**2
     +    *(efn(1,k))**2+(pnu(1,k)-Y(5,k,1)*Y(2,k,1)*fi(4,k)/
     +    Y(1,k,1))**2*(efn(4,k))**2+(fi(4,k)-Y(5,k,1)*Y(1,k,1)*
     +    pnu(1,k)/Y(2,k,1))**2*(efn(1,k))**2+(fi(1,k)-Y(5,k,1)*
     +    Y(1,k,1)*pnu(4,k)/Y(2,k,1))**2*(efn(4,k))**2)

          sum6=0
          do 65 l1=1,4
            sum6=sum6+(derf(4,1,l1))**2*(efn(l1,k))**2+
     +      (derm(4,1,l1))**2*(efn(l1,k))**2
   65     continue
 
          devY(6,k,1)=sqrt(sum6)
          
          dq12=(1/Y(8,k,1))*(di(1,2,k)-di(3,4,k))
          dq34=(-1/Y(8,k,1))*(di(1,2,k)-di(3,4,k))
          dq13=(1/Y(8,k,1))*(di(1,3,k)+di(2,4,k))
          dq24=(1/Y(8,k,1))*(di(1,3,k)+di(2,4,k))       
         
      
          sum8=0
          do 66 l1=1,4
            sum8=sum8+(dq12*derf(1,2,l1)+dq34*derf(3,4,l1)+ 
     +      dq13*derf(1,3,l1)+dq24*derf(2,4,l1))**2*(efn(l1,k))**2
     +      +(dq12*derm(1,2,l1)+dq34*derm(3,4,l1)+dq13*derm(1,3,l1)+
     +      dq24*derm(2,4,l1))**2*(efn(l1,k))**2
   66     continue
      
          devY(8,k,1)=sqrt(sum8)
           
          d7Q=-(di(4,1,k)-di(2,3,k))/(Y(8,k,1))**2
          d7d41=1/Y(8,k,1)
          d7d23=-1/Y(8,k,1)
          sum7=0

          do 67 l1=1,4

            sum7=sum7+(d7q*(dq12*derf(1,2,l1)+dq34*derf(3,4,l1)+
     +      dq13*derf(1,3,l1)+dq24*derf(2,4,l1))+d7d41*derf(4,1,l1)+
     +      d7d23*derf(2,3,l1))**2*(efn(l1,k))**2
     +      +(d7q*(dq12*derm(1,2,l1)+dq34*derm(3,4,l1)+
     +      dq13*derm(1,3,l1)+dq24*derm(2,4,l1))+d7d41*derm(4,1,l1)+
     +      d7d23*derm(2,3,l1))**2*(efn(l1,k))**2

   67     continue

          devY(7,k,1)=sqrt(sum7)

C App. res, strikes, distortion parameters...  
  
          Yst(1,k,1)=(0.2/F(k))*(Y(1,k,1)**2+Y(2,k,1)**2)
          Yst(2,k,1)=atan2(Y(2,k,1),Y(1,k,1))
          Yst(3,k,1)=0.5*atan2(-fi(3,k),fi(2,k))+rot(k)*pi/180

          if (pnu(2,k).ne.0) then
             Yst(4,k,1)=0.5*atan2(-pnu(3,k),pnu(2,k))+rot(k)*pi/180
          else
             Yst(4,k,1)=0.5*asin(Sign(1.00d+000,pnu(3,k)))+rot(k)*pi/180
          end if
       
          Yst(5,k,1)=0.5*atan((di(1,2,k)-di(3,4,k))/
     +    (di(1,3,k)+di(2,4,k)))+rot(k)*pi/180
      
C          Yst(5,k,1)=0.5*atan2((di(1,2,k)-di(3,4,k)),
C     +    (di(1,3,k)+di(2,4,k)))+rot(k)*pi/180
    
  
          rt=rot(k)*pi/180
   
          Am11=(Zr(1,2,k)+Zr(2,1,k))*sin(Yst(5,k,1)-rt)*
     +    cos(Yst(5,k,1)-rt)+Zr(1,1,k)*(cos(Yst(5,k,1)-rt))**2
     +    +Zr(2,2,k)*(sin(Yst(5,k,1)-rt))**2

          Am12=(-Zr(1,1,k)+Zr(2,2,k))*sin(Yst(5,k,1)-rt)*
     +    cos(Yst(5,k,1)-rt)+Zr(1,2,k)*(cos(Yst(5,k,1)-rt))**2
     +    -Zr(2,1,k)*(sin(Yst(5,k,1)-rt))**2

          Am21=(-Zr(1,1,k)+Zr(2,2,k))*sin(Yst(5,k,1)-rt)*
     +    cos(Yst(5,k,1)-rt)-Zr(1,2,k)*(sin(Yst(5,k,1)-rt))**2
     +    +Zr(2,1,k)*(cos(Yst(5,k,1)-rt))**2

          Am22=-(Zr(1,2,k)+Zr(2,1,k))*sin(Yst(5,k,1)-rt)*
     +    cos(Yst(5,k,1)-rt)+Zr(1,1,k)*(sin(Yst(5,k,1)-rt))**2
     +    +Zr(2,2,k)*(cos(Yst(5,k,1)-rt))**2
    
          Yst(6,k,1)=(atan(Am22/Am12))
          Yst(7,k,1)=(atan(-Am11/Am21))
          
          Yst(8,k,1)=(Yst(6,k,1)+Yst(7,k,1))/2
          Yst(9,k,1)=(Yst(6,k,1)-Yst(7,k,1))/2

          do j=1,10
            if (yst(j,k,1).eq.0) then
               yst(j,k,1)=1e-6
            end if
          end do


          devSt(1,k,1)=(0.2*2/F(i))*
     +    sqrt(Y(1,k,1)**2*devY(1,k,1)**2+Y(2,k,1)**2*devY(2,k,1)**2)
          
          devSt(2,k,1)=(1/(1+(Y(2,k,1)/Y(1,k,1))**2))*
     +    sqrt((-Y(2,k,1)/Y(1,k,1)**2)**2*devY(1,k,1)**2+
     +    (1/Y(1,k,1))**2*devY(2,k,1)**2)

          devSt(3,k,1)=0.5*(1/(1+(-fi(3,k)/fi(2,k))**2))*
     +    sqrt((fi(3,k)/fi(2,k)**2)**2*efn(2,k)**2+
     +    (-1/fi(2,k))**2*efn(3,k)**2)

          if (pnu(2,k).ne.0) then
             devSt(4,k,1)=0.5*(1/(1+(-pnu(3,k)/pnu(2,k))**2))*
     +       sqrt((pnu(3,k)/pnu(2,k)**2)**2*efn(2,k)**2+
     +       (-1/pnu(2,k))**2*efn(3,k)**2)
          end if

          devSt(5,k,1)=0.5*(1/(1+((di(1,2,k)-di(3,4,k))/
     +    (di(1,3,k)+di(2,4,k)))**2))*(1/((di(1,3,k)+di(2,4,k))**2))*
     +    sqrt((derf(1,2,1)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derf(1,3,1))**2*efn(1,k)**2+
     +    (derf(1,2,2)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derf(2,4,2))**2*efn(2,k)**2+
     +    (-derf(3,4,3)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derf(1,3,3))**2*efn(3,k)**2+
     +    (-derf(3,4,4)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derf(2,4,4))**2*efn(4,k)**2+
     +    (derm(1,2,1)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derm(1,3,1))**2*efn(1,k)**2+
     +    (derm(1,2,2)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derm(2,4,2))**2*efn(2,k)**2+
     +    (-derm(3,4,3)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derm(1,3,3))**2*efn(3,k)**2+
     +    (-derm(3,4,4)*(di(1,3,k)+di(2,4,k))-
     +    (di(1,2,k)-di(3,4,k))*derm(2,4,4))**2*efn(4,k)**2)

          do 68 ii=1,10
            Y(ii,k,2)=0.
            devY(ii,k,2)=0.
            sumY(ii)=0
            Yst(ii,k,2)=0
            sumst(ii)=0
            devst(ii,k,2)=0
   68     continue  

C  Statistical computation of invariants and related parameters
C  using nrel=100 realizations, within a gaussian error margin (=error
C  of the impedances, index = 2).
   
          CALL date_and_time(realt(1),realt(2),realt(3),jdate)
      
          n=jdate(8)+jdate(7)+jdate(6)
          ll=0
          do 80 kl=1,nrel
            n=n+1
            rafi1=gasdev(n)*efn(1,k)+fi(1,k)
            rafi2=gasdev(n+1*kl)*efn(2,k)+fi(2,k)
            rafi3=gasdev(n+2*kl)*efn(3,k)+fi(3,k)
            rafi4=gasdev(n+3*kl)*efn(4,k)+fi(4,k)
            rapnu1=gasdev(n+4*kl)*efn(1,k)+pnu(1,k)
            rapnu2=gasdev(n+5*kl)*efn(2,k)+pnu(2,k)
            rapnu3=gasdev(n+6*kl)*efn(3,k)+pnu(3,k)
            rapnu4=gasdev(n+7*kl)*efn(4,k)+pnu(4,k)      

            ry(1,kl)=Yinv12(rafi4,rafi1)
            ry(2,kl)=Yinv12(rapnu4,rapnu1)
            ry(3,kl)=Yinv34(rafi2,rafi3,rY(1,kl))
            ry(4,kl)=Yinv34(rapnu2,rapnu3,rY(2,kl))
            ry(5,kl)=(rafi4*rapnu1+rafi1*rapnu4)/(rY(1,kl)*rY(2,kl))
            rY(6,kl)=(rafi4*rapnu1-rafi1*rapnu4)/(rY(1,kl)*rY(2,kl))
            rY(8,kl)=(((rafi1*rapnu2-rafi2*rapnu1-rafi3*rapnu4+rafi4*
     +      rapnu3)**2+(rafi1*rapnu3-rafi3*rapnu1+rafi2*rapnu4-rafi4*
     +      rapnu2)**2)**0.5)/(rY(1,kl)*rY(2,kl))
            rY(7,kl)=(rafi4*rapnu1-rafi1*rapnu4-rafi2*rapnu3+rafi3*
     +      rapnu2)/(rY(1,kl)*rY(2,kl)*ry(8,kl))

C ro= st1, ph=st2, st2da=st3, st2db=st4, strtw=stdist=st5,tw1=st6,tw2=st7,
C tw=st8, sh=st9

C st3, st4 and st5: orientation of each site can be added!
C twists and shears: no rotation of the site is considered 

            rstr(1,kl)=(0.2/F(k))*(rY(1,kl)**2+rY(2,kl)**2)
            rstr(2,kl)=atan2(rY(2,kl),rY(1,kl))
            rstr(3,kl)=0.5*atan2(-(rafi3),rafi2)+rot(k)*pi/180
            rstr(4,kl)=0.5*atan2(-(rapnu3),rapnu2)+rot(k)*pi/180
            rstr(5,kl)=0.5*atan2(((rafi1*rapnu2-rafi2*rapnu1)-
     +      (rafi3*rapnu4-rafi4*rapnu3)),
     +      ((rafi1*rapnu3-rafi3*rapnu1)+(rafi2*rapnu4-
     +      rafi4*rapnu2)))+rot(k)*pi/180
            
            rt=rot(k)*180/pi
            Am11=(Rafi2)*sin(rstr(5,kl)-rt)*cos(rstr(5,kl)-rt)+
     +      (rafi1+rafi3)*(cos(rstr(5,kl)-rt))**2+(Rafi1-Rafi3)*
     +      (sin(rstr(5,kl)-rt))**2
            
            Am12=(-rafi3)*sin(rstr(5,kl)-rt)*cos(rstr(5,kl)-rt)+
     +      (rafi2+rafi4)*(cos(rstr(5,kl)-rt))**2-(rafi2-rafi4)*  
     +      (sin(rstr(5,kl)-rt))**2
            
            Am21=(-rafi3)*sin(rstr(5,kl)-rt)*cos(rstr(5,kl)-rt)-
     +      (rafi2+rafi4)*(sin(rstr(5,kl)-rt))**2+(rafi2-rafi4)*
     +      (cos(rstr(5,kl)-rt))**2      
            
            Am22=(-rafi2)*sin(rstr(5,kl)-rt)*cos(rstr(5,kl)-rt)+
     +      (rafi1+rafi3)*(sin(rstr(5,kl)-rt))**2+(rafi1-rafi3)*
     +      (cos(rstr(5,kl)-rt))**2
            
            rstr(6,kl)=-atan(Am22/Am12)
            rstr(7,kl)=-atan(-Am11/Am21)
            
            rstr(8,kl)=(rstr(6,kl)+rstr(7,kl))/2
            rstr(9,kl)=(rstr(6,kl)-rstr(7,kl))/2
            
            do 81 ii=1,9
              sumY(ii)=sumY(ii)+rY(ii,kl)
              sumst(ii)=sumst(ii)+rstr(ii,kl)
   81       continue
   80     continue
      
          do 82 ii=1,9
            Y(ii,k,2)=sumY(ii)/nrel
            dify(ii)=0
            Yst(ii,k,2)=sumst(ii)/nrel
            difst(ii)=0
   82     continue
         
            do 83 la=1,nrel
            do 84 ii=1,10
              dify(ii)=dify(ii)+(Y(ii,k,2)-ry(ii,la))**2
              difst(ii)=difst(ii)+(Yst(ii,k,2)-rstr(ii,la))**2
   84       continue
   83     continue
          
          do 85 ii=1,9
            devY(ii,k,2)=sqrt(dify(ii)/(nrel-1))
            devst(ii,k,2)=sqrt(difst(ii)/(nrel-1))
   85     continue

          do 86 ii=1,8
            Y(ii,k,1)=abs(Y(ii,k,1))
            Y(ii,k,2)=abs(Y(ii,k,2))
   86     continue
   60   continue

C     Using true or statistic values:

C     true values: 
        li=1

C     statistic values 
C     li =2

C     Writing output files in frequency decreasing order:      
        if (f(1).lt.f(2)) then
           finc=-1
          else
             finc=1
        end if 

        do 90 ka=1,nfreq
          if (finc.eq.-1) then
             k=nfreq+1-ka
          else
             k=ka
          end if

          ig=0      
          
          do 94 i=1,8
            bias(i)=y(i,k,1)-y(i,k,2)
   94     continue
      
          do 941 i=3,6

            if (abs(bias(i)).gt.devY(i,k,2)) then
               if ((y(i,k,1).lt.th).and.
     +         (y(i,k,2).lt.th).and.
     +         (abs(bias(i)).lt.th)) then
                  ig=ig+0
               else
                  ig=ig+1
               end if
            else
               ig=ig+0
            end if
  941     continue
 
          do 944 i=2,9
            yst(i,k,1)=yst(i,k,1)*180./pi
            yst(i,k,2)=yst(i,k,2)*180./pi
            devst(i,k,1)=devst(i,k,1)*180/pi
            devst(i,k,2)=devst(i,k,2)*180/pi
  944     continue

C     first quadrant of strike angles
          do 945 i=3,5
            yst(i,k,1)=fq(yst(i,k,1))
            yst(i,k,2)=fq(yst(i,k,2))
  945     continue 
  
          igst=0
          do 96 i=1,9
            biasst(i)=yst(i,k,1)-yst(i,k,2)
            if (abs(biasst(i)).gt.devst(i,k,2)) then
               igst=ig+1.
            else
               igst=ig+0
            end if
   96     continue 

          do 955 j=3,7
            yc(j,k)=34000
  955     continue

          if (ig.eq.0) then
             if (((y(8,k,li).le.thq).or.(y(7,k,li).gt.1.)).or.
     +       ((y(8,k,li).le.thq).and.(y(7,k,li).gt.1.))) then
                yc(7,k)=1000
             else
                if (y(7,k,li).gt.th) then
                   yc(7,k)=1
                else 
                   yc(7,k)=0
                end if
             end if

             do 91 j=3,6
               tin=0
               if ((y(j,k,li)+devY(j,k,li).gt.th).and.
     +         (y(j,k,li)-devY(j,k,li).gt.th)) then
                  yc(j,k)=1
                  tin=1
               end if
               if ((y(j,k,li)+devY(j,k,li).le.th).and.
     +         (abs(y(j,k,li)-devY(j,k,li)).le.th)) then
                  yc(j,k)=0
                  tin=1
               end if
               if (tin.eq.0) then
                  if ((y(j,k,li)+devY(j,k,li).gt.th)
     +            .and.(y(j,k,li).gt.th).and.
     +            (abs(y(j,k,li)-devy(j,k,li)).le.th)) then
                     yc(j,k)=1
                  end if
                  if ((y(j,k,li)+devy(j,k,li).gt.th).and.(y(j,k,li)+
     +            devy(j,k,li).le.1).and.(y(j,k,li).le.th).and.
     +            (abs(y(j,k,li)-devy(j,k,li)).le.th))then
                     yc(j,k)=1
                  end if
               end if
   91        continue
          end if
     
 2003     t=0
          do 92 ja=3,7
            if (yc(ja,k).ne.34000) then
               t=t+1
               nb=nb+1
            end if
   92     continue
          if (t.eq.5) ta=ta+1

C Dimensionality classification for each frequency:
          icas=0
          if ((yc(3,k).ne.34000).and.(yc(4,k).ne.34000).and.
     +    (yc(5,k).ne.34000).and.(yc(6,k).ne.34000).and.
     +    (yc(7,k).ne.34000)) then

C 1D    
             if((yc(3,k)+yc(4,k)+yc(5,k)+yc(6,k).eq.0).and.   
     +       ((yc(7,k).eq.1000).or.(yc(7,k).eq.0))) then
               icas=1
             end if

C distortion over 1D o 2D, o 2D:
             if (((yc(3,k)+yc(4,k)).ge.1).and.(yc(5,k)+yc(6,k).eq.0)
     +       .and.((yc(7,k).eq.1000).or.(yc(7,k).eq.0))) then
                if ((abs(0.5*(Zr(1,2,k)-Zr(2,1,k))).lt.th).and.
     +          (abs(0.5*(Zi(1,2,k)-Zi(2,1,k))).lt.th)) then
                   icas=6
                else
                   if ((abs(yst(3,k,2)-yst(4,k,2)).lt.10).or.
     +             (abs(yst(3,k,2)-yst(4,k,2)).gt.(80))) then
                      icas=2
                   else
                      icas=4
                   end if
                end if
             end if  

C distortion twist

             if (((yc(3,k)+yc(4,k)).ge.1).and.(yc(5,k).eq.1).and.
     +       (yc(6,k).eq.0).and.(yc(7,k).eq.0)) then
                icas=3
             end if
             
C distortion
             if (yc(6,k).eq.1) then
                if (yc(7,k).eq.0) then
                   icas=4
                end if
             end if

C 3D
             if (yc(7,k).eq.1) then
                icas=5
             end if
 
C 3D/1D2D      
             if (((yc(3,k)+yc(4,k)).ge.1).and.(yc(5,k).eq.1)
     +       .and.(yc(6,k).eq.0).and.(yc(7,k).eq.1000)) then
                icas=7
             end if
          end if


C  Converting dimensionality cases (integer) into character strings:



          CALL WriteFile3(mp,k,plongi(l),plati(l),sidsite(l),
     +                   f(k),Y,devY)

          CALL WriteFile4(mp,k,sidsite(l),plongi(l),
     +                   plati(l),f(k),icas,yst,devst)
 
          CALL WriteFile7(mp,k,li,sidsite(l),plongi(l),
     +                   plati(l),f(k),rot(k),icas,y,devy,yst,devst)
          
            
            if ((statf.eq.'y').or.(statf.eq.'Y')) then
             CALL WriteFile8(mp,k,sidsite(l),plongi(l),
     +                      plati(l),f(k),y,devy,yst,devst,bias,biasst)
          end if
                        
            CALL WriteFile9(mp,k,sidsite(l),plongi(l),
     +                   plati(l),f(k),icas,yst,devst,namedim)

          CALL CompWriteFile10(mp,sidsite(l),k,f(k),
     +                   zr,zi,fi,pnu,y,yst,icas,namedim)
         
          if (icas.lt.12) then
             nicas(icas)=nicas(icas)+1
          end if
          if (icas.eq.2) then
             nn=nn+1
             dst34=dst34+abs(yst(3,k,1)-yst(4,k,1))
          end if

C     Save variables needed for band averaging:

          if ((avg.eq.'y').or.(avg.eq.'Y')) then
      
             af(ka,l)=f(k)
               iacas(ka,l)=icas
             arot(ka,l)=rot(k)
             
             do i=1,9
               ayst(i,ka,l)=yst(i,k,1)
               adevst(i,ka,l)=devst(i,k,2)
             end do
          end if
   
C reset dimensionalities:
          icas=0


   90   continue
 
        goto 204

  201   write(6,202) file2(1:lsidsite(l)+4)
  
  202 format(' File ',a,' not found')

        write(6,*) ' Press any key to continue running the program or
     +STOP to finish'
        read(5,'(a)') info
          write(6,*)
          if ((info.eq.'STOP').or.(info.eq.'stop')) then
             goto 700
          end if

  101   write(6,*) 'WARNING: In file, ',file11,'the 
     +number of sites is higher than the actual list of edi files'
        write(6,*) 'Press any key to continue or STOP to finish'
        read(5,'(a)') info
          write(6,*)
          if ((info.eq.'STOP').or.(info.eq.'stop')) then
             goto 700
          end if

        lmiss=lmiss+1
         
        goto 20

  204   close(9)
   20 continue
   
  203 close(3)
      close(11)

      close(4)
      close(7) 
      close(8)
      close(10)
  
      nsit=nsites-lmiss



      if ((nsit).ne.0) then

       write(6,*) ' ***************************************************'
       write(6,*)
       write(6,*) ' Dimensionality analysis for each site and frequency'
       write(6,*) ' finished'
       write(6,*)


       write(6,111) nsit
  111 format('  Dimensionality summary from the ',i3,' EDI files:')
       write(6,*)
       write(6,*) ' Total periods=',lc
       write(6,*) ' Dimensionality cases:'
       write(6,*) ' Undetermined cases=   ', nicas(0)
       write(6,*) ' 1D cases=             ', nicas(1)
       write(6,*) ' 2D cases=             ', nicas(2)
       write(6,*) ' 3D/2Dtwist cases=     ', nicas(3)
       write(6,*) ' 3D/2D cases=          ', nicas(4)
       write(6,*) ' 3D cases=             ', nicas(5)
       write(6,*) ' 3D/2Ddiag cases=      ', nicas(6)
       write(6,*) ' 3D/1D2D cases=        ', nicas(7)
       write(6,*)
       write(6,*) ' ***************************************************'
  
      else
      
       write (6,*) 'No parameters and/or dimensionality analysis could '
       write (6,*) 'be performed from any of the sites in the list.'
       write (6,*)
       write (6,*) 'Check list with files, edi files formats or file 
     + param.cfg'
       write (6,*) 'Program will finish'

       goto 700
      end if

C------------------------------------------------------------------------------

C  Grouping in bands:
      
      if ((avg.eq.'y').or.(avg.eq.'Y')) then

         dpm=log10(pm)
         dpma=log10(pma)
         
         pinc=1./bpd
         
         ng=(dpma-dpm)*bpd
         
         write(6,*)
         write(6,*) 'Averaging results in period bands '
         write(6,141) pm
  141 format(' Min per =  ', e10.4)
           write(6,142) pma
  142 format (' Max per =  ', e10.4)
         write(6,*) ng, ' bands'      
         write(6,*) 

         do kl=1,ng

               perm=10**(dpm+(kl-1)*pinc)
               permA=10**(dpm+kl*pinc)
               
               write(6,601) kl,perm, perma
  601 format('Band',i2,':  ',g10.4,' s - ',g10.4,' s')
               write(6,*)
            
         end do


         count1=0

         do 500 nc=1,ng
  
           count2=0

           do 510 ii=1,nsit
             
             if (nc.eq.1) then
             open(unit=9,file=file9(ii),status='old')
             
                do il=1,6+nf(ii)
                  read(9,'(a)') info
                end do
                
                CallHeader9bis
              end if


             i1=0
             i2=0
             i1old=0
             i2old=0
             
             ifirst=0
             
             lsi=ii
             
             do 490 kl=1,ng

               ia1=0      
               
               permin=10**(dpm+(kl-1)*pinc)
               permax=10**(dpm+kl*pinc)
               
               do 540 i=1,Nf(lsi)

                 per(i)=1/af(i,lsi)
                 
                 if ((per(i).gt.permin).and.
     +           (per(i).le.permax)) then
                    ifirst=ifirst+1
                    ia1=ia1+1
                    if (ia1.eq.1) then
                       if (ifirst.eq.1) then
                          i1=i
                          i2=i
                       else
                          i1=i2+1
                          i2=i2+1
                       endif
                    else
                       i2=i2+1
                    end if
                 end if
  540          continue
      
               do l=1,7
                 ict(l)=0
               end do
               maxi=0                
               pin=0
               ipin=0

               IF ((i2.ge.i1).and.(i2.ne.0)
     +         .and.(i2.ne.i2old).and.(i1.ne.i1old)) THEN

                  do 600 i=i1,i2
                    if (iacas(i,lsi).ne.0) then
                       pin=pin+1
                       ipin=ipin+1
                       if (iacas(i,lsi).eq.3) then
                          iacas(i,lsi)=4
                       end if
                       do j=1,7
                         if (iacas(i,lsi).eq.j) then
                            ict(j)=ict(j)+1 
                         end if
                       end do
                    end if
  600             continue
                  
c  Mode:

                  maxx=ict(1)
                  do la=2,7
                    maxx=max(maxx,ict(la))
                  end do
                  
                  strike=0
                  shear=0
                  twi=0
                  j2=0
                  errstr=0
                  errtwi=0
                  errsh=0

                  scale=0
                  
                  if ((maxx.eq.0).or.(pin.le.1)) then
                     goto 2001
                  else
      
                     do la=7,1,-1
                       if (ict(la).eq.maxx) then
                          maxi=la
                       else
                          do lb=la+1,7
                            if ((ict(lb).eq.maxx).and.
     +                      (ict(la).eq.maxx)) then
                               maxi=min(la,lb)
                            end if 
                          end do
                       end if 
                     end do
      
                     if (maxi.eq.5) then
      
                        do lc=7,1,-1
                          if ((lc.ne.5).and.(ict(lc).ne.0)) then
                             if (ict(lc).ge.pin/2) then
                                maxi=lc 
                             end if
                          end if
                        end do
                     end if

                  end if

                 if (maxi.eq.2) then
                    do 610 i=i1,i2
                      if (iacas(i,lsi).eq.maxi) then
                         j2=j2+1

C Angles average:
                         if(abs(ayst(3,i,lsi)-ayst(4,i,lsi)).gt.80) then
                            if (ayst(4,i,2).lt.ayst(3,i,lsi)) then
                               ayst(4,i,lsi)=90-ayst(4,i,lsi)
                            else
                               ayst(3,i,lsi)=90-ayst(3,i,lsi)
                            end if
                         end if
                       strike=strike+(ayst(3,i,lsi)+ayst(4,i,lsi))/2
                       errstr=errstr+((adevst(3,i,lsi))**2+
     +                 (adevst(4,i,lsi))**2)/4
                      end if
  610               continue
                    strike=strike/j2
                    errstr=sqrt(errstr)/j2
                    scale=1+1/errstr
                 end if
      
      
                 if (maxi.eq.4) then
                    do 630 i=i1,i2
                      if (iacas(i,lsi).eq.maxi) then
                         j2=j2+1
                         if (i.ne.i1) then
                            if (((aYst(5,i-1,lsi).ge.80).and.
     +                      (aYst(5,i,lsi).le.10)).or.
     +                      ((aYst(5,i-1,lsi).le.10).and.
     +                      (aYst(5,i,lsi).ge.80))) then
                               ayst(5,i,lsi)=90-ayst(5,i,lsi)
                            end if
                         end if
         
                         strike=strike+aYst(5,i,lsi)
                         errstr=errstr+(adevst(5,i,lsi))**2
                         twi=twi+aYst(8,i,lsi)
                         errtwi=errtwi+(adevst(8,i,lsi))**2
                         shear=shear+aYst(9,i,lsi)
                         errsh=errsh+(adevst(9,i,lsi))**2
                      end if
  630               continue

                    strike=strike/j2
                    shear=shear/j2
                    twi=twi/j2
                    errstr=sqrt(errstr)/j2
                    errsh=sqrt(errsh)/j2
                    errtwi=sqrt(errtwi)/j2
                    scale=1+1/errstr

                 end if
     
 2001            if ((maxi.ge.2).and.(maxi.le.4)) then
                    nst=1
                 else
                    nst=0
                 end if


                 if (nc.eq.1) then

                    count1=count1+1

                  if (count1.eq.1) then
                     CALL Outputheader14(output,nout,th)
                    end if

                    CALL WriteFile14(mp,sidsite(lsi),
     +                 plongi(lsi),plati(lsi),kl,per,ipin,i1,i2,maxi,
     +                 strike,errstr,twi,errtwi,shear,errsh,
     +                 nst,scale)
      

                    CALL writefile9bis(mp,sidsite(lsi),
     +                 plongi(lsi),plati(lsi),kl,per,ipin,
     +                 i1,i2,namedim,maxi,strike,errstr,twi,errtwi,
     +                 shear,errsh,nst,scale)
                 end if

                 if (kl.eq.nc) then

                    count2=count2+1
     
                  if (count2.eq.1) then
                     
                   CALL Header15(output,nout,nc,th,permin,permax)
                  end if

                    write(15,151) sidsite(lsi),plongi(lsi),plati(lsi),
     +              ipin,maxi,strike,errstr,twi,errtwi,shear,errsh,nst,
     +              scale,strike*(-1)
  151 format(a,2(g14.6,2x),i4,2x,i3,1x,6(f9.4),2x,i2,2(f9.4))
 
                 end if

                 ELSE
 
                 if ((nc.eq.1)) then
                    write(9,901) sidsite(ii),plongi(ii),plati(ii),kl
                    
  901 format(a20,2(g14.6,2x),
     +1x,i3,2x,'NO PERIODS AVAILABLE FOR THIS BAND') 
                    
                    write(6,602) kl,sidsite(ii)
  602 format('NO PERIODS AVAILABLE FOR BAND ',i2,' FOR SITE ',
     +a)
                    write(6,*) 

                 end if 
 
               END IF

            i1old=i1
            i2old=i2

  490       continue

            close(9)

  510     continue
   
          close(15)
  
  500   continue
  
        close(14)

      end if
  
C------------------------------------

  699 write(6,*) 'Program finished, bye.'
      

  700 stop
  
      end