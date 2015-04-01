        subroutine getsp1(nch,sampr,scale,nfil,iftype,afparam,iounit
     1   ,decl,stcor,orient,chid,cda,cdb,lclkd)
        
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

        include 'fcor.inc'
        
        real scale(*),orient(2,*),stcor(2),decl
        integer iftype(nfilmax,*),nfil(*)

        character*80 afparam(nfilmax,*),cfsp
        character*40 ctemp
        character*2 ftype               
        character*1 chid(*)
        character*20 cname
        logical lclkd

c  station name
c        read(iounit,*)
        read(iounit,'(a20)',err=100) cname
        print*,'cname',cname
c   position of station; latitude and longitude; fractions of degrees
c   are assumed expressed as decimals (i.e. 45:30:00 = 45.5000)
        read(iounit,*,err=100) stcor      
          print*,'stcor',stcor
c   magnetic declination of station
        read(iounit,*,err=100) decl
          print*,'decl',decl
c    number of data channels
        read(iounit,*,err=100) nch  
         if(nch.gt.nchmx) then
            write(0,*) 'ERROR: nchmx = ',nchmx,' is not big enough'
            write(0,*) '   nch = ',nch,' ; increase nchmx and recompile'
            stop
         endif
           print*,'nch',nch
c     sampling rate; assumed in seconds
        read(iounit,*,err=100) sampr
           print*,'sampr',sampr
c   clock drift corrections; assume clock time = tc; actual time = t
c     these are related by tc = t + (cda + cdb*t); assumes t in seconds
        read(iounit,*) cda,cdb
           print*,'cda,cdb',cda,cdb
        lclkd = ( abs(cdb) .gt. 1.0e-20)
           
c      get filter parameters
           
           do 10 i = 1,nch

           read(iounit,'(a1)',err=100) chid(i) 
               print*,'chid(i) ',i,chid(i)

           if( (chid(i).eq.'e') .or. (chid(i).eq.'E') ) then
c      electric field channel:
c        read in electrode line length, angle (deg.  e of geomangnetic north of line
c          from neg. to pos. electrode), amplifier gain
              read(iounit,*,err=100) r,orient(1,i),orient(2,i),ampg

c       read in count conversion ( = nt/count);
c       number of filters to correct for for electric field channels
              read(iounit,*,err=100) scale(i),nfil(i)
              scale(i) = scale(i)/(ampg*r)

           else

c     magnetic field channel (or something else ...):
c       assume geomagnetic orientation of measurements
c       read in count conversion ( = nt/count);
c       number of filters to correct for for magnetic field channels
              read(iounit,*) orient(1,i),orient(2,i)
                 print*,orient(1,i),orient(2,i)
              read(iounit,*) scale(i),nfil(i)
                 print*,'scale,nfil',scale(i),nfil(i)
           end if

c>>>>>>>>>> get analogue filter/system response info
              do 8 j = 1,nfil(i)
              read(iounit,'(a2)',err=100) ftype
              read(iounit,'(a80)',err=100) afparam(j,i)
              
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
              else 
                 iftype(j,i) = 9
                 write(6,*) 'filter type',ftype,' unknown'
              end if
            
8             continue

10         continue

        close(iounit)
        return

100     continue    

c       error trying to read system parameter file (probably doesn't
c       exist)

           print*,'error reading system parameter file'
           stop

        end
