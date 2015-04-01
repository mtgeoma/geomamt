        subroutine getsp(nch,sampr,scale,nfil,iftype,afparam 
     1   ,decl,stcor,orient,chid,cda,cdb,lclkd)
        
      include 'iounits.inc'

c       reads info from file to set up filter corrections, conversion
c       of counts to physical units, conversion of elctrode coordinates
c       to an orthogonal (geomagnetic) coordinate system
c
c        inputs required:
c              nch = # of channels
c              dr =  sampling rate (needed only if no system parameter
c                                 file is provided.... for gds data only
c                               the routine will assume a two pole low pass
c                               filter with a time constant equal to the 
c                               nyquist period.)
c               stname =   character string giving station name; routine
c                         trys to open a system parameter file named spstname

c>>>>>>>>>>>> 28 Feb, 1991: afparam is changed to an aray of
c       character*80 variables to allow more general filter/response
c         parameters (e.g., calibration files)

cf 7.6.95 Magson Magnetometer 10 sec filter added

cnew    Add two filter identifiers for coefficients of analog filters
cnew    provided in form of simple tables, one as real and imaginary part
cnew    and one as amplitude and phase, both as functions of period.
cnew    Within the table files all header lines should be marked with one of the
cnew    standard comment characters such as #, %, $, &, *
cnew    AP is the identifier for the amplitude/phase table type
cnew    RI is the identifier for the rela/imaginary table type

        include 'params2.inc'
        
        real scale(*),orient(2,*),stcor(2),decl
        integer iftype(nfilmax,*),nfil(*)

        character*80 afparam(nfilmax,*)
        character*2 ftype               
        character*6 chid(*)
        character*20 cname
        character*20 ctemp
        logical lclkd

        open(unit =sp_unit, file = cfsp,status='old')
c  station name
        read(sp_unit,'(a20)',err=100) cname
c   position of station; latitude and longitude; fractions of degrees
c   are assumed expressed as decimals (i.e. 45:30:00 = 45.5000)
        read(sp_unit,*,err=100) stcor      
c   magnetic declination of station
        read(sp_unit,*,err=100) decl
c    number of data channels
        read(sp_unit,*,err=100) nch  
c     sampling rate; assumed in seconds
        read(sp_unit,*,err=100) sampr
c   clock drift corrections; assume clock time = tc; actual time = t
c     these are related by tc = t + (cda + cdb*t); assumes t in seconds
        read(sp_unit,*) cda,cdb
        lclkd = ( abs(cdb) .gt. 1.0e-20)

        write (6,22) cname,nch, 1./sampr, stcor, decl 
 22     format(/,'Station:',a9,3x,'# of channels:',i3,3x,' Dt:',f8.2,
     &    'Hz',/,'Lat./Long.:',2f7.2,3x,' Declination:',f6.2,/)

c      get filter parameters
        write(6,23)
 23     format('    CHID    ORIENT.   TILT    SCALE          FILTERS')
           do 10 i = 1,nch

           read(sp_unit,'(a6)',err=100) chid(i) 

           if( (chid(i)(1:1).eq.'e') .or. (chid(i)(1:1).eq.'E') ) then
c      electric field channel:
c        read in electrode line length, angle (deg.  e of geomangnetic north of line
c          from neg. to pos. electrode), amplifier gain
              read(sp_unit,*,err=100) r,orient(1,i),orient(2,i),ampg

c       read in count conversion ( = nt/count);
c       number of filters to correct for for electric field channels
              read(sp_unit,*,err=100) scale(i),nfil(i)
              scale(i) = scale(i)/(ampg*r)

           else

c     magnetic field channel (or something else ...):
c       assume geomagnetic orientation of measurements
c       read in count conversion ( = nt/count);
c       number of filters to correct for for magnetic field channels
              read(sp_unit,*) orient(1,i),orient(2,i)
              read(sp_unit,*) scale(i),nfil(i)
           end if


c>>>>>>>>>> get analogue filter/system response info
           ctemp = ' '
           do 8 j = 1,nfil(i)
              read(sp_unit,'(a2)',err=100) ftype
              ctemp = ctemp(1:iclong(ctemp,20))//ftype//' '
              read(sp_unit,'(a80)',err=100) afparam(j,i)
              
              if((ftype.eq.'l1').or.(ftype.eq.'L1')) then
c ..................one pole low pass filter
c                    2 real parameters:
                 iftype(j,i) = 5
              else if((ftype.eq.'l2').or.(ftype.eq.'L2')) then
c ..................two pole low pass filter
c                    3 real parameters:
                 iftype(j,i) = 4
              else if((ftype.eq.'h1').or.(ftype.eq.'H1')) then
c ..................one pole hi pass filter
c                    2 real parameters:
                 iftype(j,i) = 2
              else if((ftype.eq.'h2').or.(ftype.eq.'H2')) then
c ..................two pole hi pass filter
c                    3 real parameters:
                 iftype(j,i) = 3
              else if((ftype.eq.'bc').or.(ftype.eq.'BC')) then
c .................."box car" average
c                    3 real parameters:
                 iftype(j,i) = 6
              else if((ftype.eq.'te').or.(ftype.eq.'TE')) then
c ...................EMI electric field calibration table
c                         Three parameters:
c                         file name, mode of operation, gain
                 iftype(j,i) = 7
              else if((ftype.eq.'tb').or.(ftype.eq.'TB')) then
c ...................EMI mag field calibration table
c                         One parameter:
c                         file name
                 iftype(j,i) = 8

C-PDAS  
C-PDAS  31/8/96 these are the old analytical transfer functions
C-PDAS  new ones, based on the true calibration files have names
C-PDAS  PA/pa and 05 for PAB and MFS05, respect., see below
 
              else if ((ftype.eq.'pd').or.(ftype.eq.'PD')) then
c....................PDAS 200 Hz low pass filter
c                    all parameters in FUNCTION AFCOR
                iftype(j,i) = 11
              else if ((ftype.eq.'mx').or.(ftype.eq.'MX')) then
c....................METRONIX GMS05 coil
c                    all parameters in FUNCTION AFCOR
              iftype(j,i) = 12

C-PDAS end
C  LMT 
              else if ((ftype.eq.'RA').or.(ftype.eq.'ra')) then
c....................10 sec low pass of RAP telluric
c                    all parameters in FUNCTION AFCOR
                iftype(j,i) = 13
              else if ((ftype.eq.'MA').or.(ftype.eq.'ma')) then
c....................MAGSON transfer function (5 sec)
c                    all parameters in FUNCTION AFCOR
                iftype(j,i) = 14
              else if ((ftype.eq.'R2').or.(ftype.eq.'r2')) then
c....................120 sec low pass of RAP telluric
c                    all parameters in FUNKTION AFCOR
                iftype(j,i) = 15
              else if ((ftype.eq.'M1').or.(ftype.eq.'m1')) then
c....................MAGSON transfer function (10 sec)
c                    all parameters in FUNCTION AFCOR
                iftype(j,i) = 16
C  END LMT

C-PDAS  new transfe functions based on calibration data files
              else if ((ftype.eq.'PA').or.(ftype.eq.'pa')) then
C....................PAB KMT preamp, data in calibration file
                iftype(j,i) = 17
              else if (ftype.eq.'05') then
C..................MFS05 transfer functions from calibration data file
                iftype(j,i) = 18
C-PDAS   new end

c LRMT Bessel filter:
             else if ((ftype.eq.'BS').or.(ftype.eq.'bs')) then
                iftype(j,i) = 19

c metronix mfs transfer function (theoretical and/or calibration files):
             else if ((ftype.eq.'MS').or.(ftype.eq.'ms')) then
                iftype(j,i) = 20

             elseif (ftype.eq.'UC') then
ccc             uncallibrated coil i omega response
                iftype(j,i) = 97
             elseif (ftype.eq.'AP') then
                iftype(j,i) = 98
             elseif (ftype.eq.'RI') then
                iftype(j,i) = 99

             else 
                 iftype(j,i) = 9
                 write(6,*) 'filter type',ftype,' unknown'
              end if
            
8             continue

           write (6,'(i2,2x,a6,5x,2f7.2,e12.5,i3,1x,a9)')i,chid(i),
     &       orient(1,i),orient(2,i),scale(i),nfil(i),
     &       ctemp(1:iclong(ctemp,20))

10         continue

        close(sp_unit)
        return

100     continue    

c       error trying to read system parameter file (probably doesn't
c       exist)

           print*,'error reading system parameter file'
           stop

        end
