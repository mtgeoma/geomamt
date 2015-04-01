c----------------------------------------------------------------------

      PROGRAM EDI2J 

c---------------------------------------------------------------------
c Versions
c  1.0       Original
c  1.1  agj  Edited resread.f for missing ZROT in IMPEDANCE EDI file
c            Corrected Z error estimation for IMPEDANCE EDI file
c            Writes in increasing period
c            Checks whether TZ in file, and only writes out if existing
c 1.2 agj Does not write out MISDAT data into J-file
c 1.3 agj Can state Spectra, Impedance or Unknown
c
c---------------------------------------------------------------------

      implicit none

c---------------------------------------------------------------------

      include 'ctrlblk.inc'

c---------------------------------------------------------------------

c...if you change MAXF here, you must also change it in edispecread and editfread !!!
      integer MAXF
      parameter ( MAXF = 200 )
      real PI
      parameter ( PI = 3.141592 )
      real MU
      parameter ( MU = 4.*PI*1.E-7 )

c---------------------------------------------------------------------

      COMPLEX Z(3,2), zdata(6,MAXF), zz
      real zvar(6,MAXF)
      REAL rho(2,2), phs(2,2), srh(2,2), sph(2,2)
      REAL vimp(3,2), delta, avgs, freq(MAXF), angle/0.0/
      REAL sdata(7,7,MAXF), info(MAXF,3), temp(7,7)
      REAL rxx(MAXF), ryy(MAXF), rxy(MAXF), ryx(MAXF)
      REAL pxx(MAXF), pyy(MAXF), pxy(MAXF), pyx(MAXF)
      REAL zxxr(MAXF), zyyr(MAXF), zxyr(MAXF), zyxr(MAXF)
      REAL zxxi(MAXF), zyyi(MAXF), zxyi(MAXF), zyxi(MAXF)
      REAL tzxr(MAXF), tzyr(MAXF), tzxi(MAXF), tzyi(MAXF)
      REAL vzxx(MAXF), vzyy(MAXF), vzxy(MAXF), vzyx(MAXF)
      REAL vtzy(MAXF), vtzx(MAXF), one/1.00/
      REAL srxx(MAXF), sryy(MAXF), srxy(MAXF), sryx(MAXF)
      REAL spxx(MAXF), spyy(MAXF), spxy(MAXF), spyx(MAXF)
      real zerr(4,MAXF), zzerr, period, rhoa, phase, zr, zi, phaerr
      real rhomax, rhomin, phamax, phamin
      real datang
      character ifile*80/'wml-03nb.edi'/, ofile*80
      character site*6, list*40, sitefile*40
      character*24 date_string
      character*80 line
      character*4 ver
      character*1 char1
      integer nfreq, i, j, k, lnblnk, leq
      integer istart, iend, iinc, itz
      integer imiss(MAXF)

      real ddmmss2deg
      real xlat, xlon, elev
      real misdat
      integer nrhomisdat, ntzmisdat
      logical SPECTRA
      logical ONEEDI
      logical EGBERT
      logical UNKNOWN
      logical AUTONAME

      real convz2r, convz2p, factor

      real cosd, sind, atand, atan2d
      external cosd, sind, atand, atan2d

c-------------------------------------------------------------------------------
c...initialize

c...factor to convert impedances from field units to S.I.
c   This is mu0 (to convert from B in nT to H in gamma) 
c   divided by 10^-3 (to convert from mV/km to V/m)
      factor = 4.*PI*1.e-4

c-------------------------------------------------------------------------------
c...set version number: This must be the same number as "version.mk"      
      ver = 'v1.2'

c-------------------------------------------------------------------------------
c...get iitialize some parameters
      do i = 1, MAXF
        imiss(i)=1
      enddo

      misdat = 0
      nrhomisdat = 0

c-------------------------------------------------------------------------------
c...get date

      call fdate(date_string)

      write(*,*)'        Welcome to edi2j '
      write(*,*)'         version = ', ver
      write(*,*)' Todays date: ', date_string
      write(*,*)

c-------------------------------------------------------------------------------
c...get print value

      iprint = 0
      call iin('Give print level (0/1/2/3)', iprint )

c-------------------------------------------------------------------------------
c...get name of edi file to read

      list = 'edi2j.lst'
      call tin('Enter EDI-file list OR single file (.edi)',list)
      if( index(list,'.edi').eq.0 .and. index(list,'.EDI').eq.0) then
        ONEEDI = .FALSE.
        open(unit=2,file=list,status='old')
      else
        ONEEDI = .TRUE.
      endif

      EGBERT = .FALSE.
c      call lin('Is this an EGBERT EDI file?', EGBERT )

      AUTONAME = .TRUE.
      call lin('Automatically name the site?', AUTONAME )

      char1 = 'U'
      call tin('Are data in Spectra(S), Impedance(I) or Unknown(U)', 
     &          char1 )
      call shiftup( char1 )
      if( char1.eq.'S' ) then
        SPECTRA = .TRUE.
        UNKNOWN = .FALSE.
      elseif( char1.eq.'I' ) then
        SPECTRA = .FALSE.
        UNKNOWN = .FALSE.
      elseif( char1.eq.'U' ) then
        SPECTRA = .FALSE.
        UNKNOWN = .TRUE.
      else
        write(*,*)'Unknown response, must be S, I or U'
        stop
      endif

      do while(.TRUE.)

        if( ONEEDI ) then
          ifile = list
        else
          read(2,'(a80)',end=999) ifile
        endif  

        if( iprint.ge.1 ) then
          write(*,*)'Opening input EDI file >', ifile(1:lnblnk(ifile))
        endif
        sitefile = ifile(1:index(ifile,'.edi')-1)
c...get rid of directories in sitefile
        do while( index(sitefile,'/').ne.0 ) 
          sitefile = sitefile(index(sitefile,'/')+1:)
        enddo
        write(*,'(2a)')'>>>>>Sitefile = ', sitefile(1:lnblnk(sitefile))

c...read edi-file - if UNKNOWN search file for keyword ">=SPECTRASECT"
        if( UNKNOWN ) then
          open(unit=4,file=ifile)
          read(4,'(a)') line
          do while ( index(line,'>=SPECTRASECT').eq.0 )
            read(4,'(a)',end=3 ) line
          enddo
          SPECTRA = .TRUE.
          if( iprint.ge.0 )write(*,*) 'Opening SPECTRA file >', 
     &                             ifile(1:lnblnk(ifile))
          goto 4
3         SPECTRA = .FALSE.
          if( iprint.ge.0 )write(*,*) 'Opening IMPED   file >', 
     &                             ifile(1:lnblnk(ifile))
4         close(unit=4)
        endif

      if( SPECTRA ) then
        call edispecread(ifile,misdat,freq,nfreq,sdata,info,site)
      else
        call editfread(ifile,misdat,freq,nfreq,zdata,zvar,info,site)
      endif

      if( .not.AUTONAME ) then
        call tin('Give site name (6 char)', site )
      endif

      if( nfreq.gt.MAXF ) then
        write(*,*)'Number of frequencies too great for this version'
        write(*,*)'MAXF =', MAXF,'    nfreq =', nfreq
        stop
      endif

      if(iprint.ge.1) then
        write(*,*)'nfreq read in =', nfreq
      endif

      datang = info(1,1)
      do i = 2, nfreq
        if( info(i,1).ne.datang ) then
          write(*,*)'Not set up for different angles at each freq'
          stop
        endif
      enddo

      if( iprint.ge.1 ) write(*,*)'Data in angle of:', datang
      angle = datang
c      call rin('Angle to rotate data to', angle)


      delta = angle - datang
c      write(*,*) 'Data in angle of         :', datang
c      write(*,*) 'Data required in angle of:', angle
c      write(*,*) 'Rotation applied         :', delta

c...open output file

      ofile = sitefile(1:lnblnk(sitefile))//'.dat'
      call shiftdown(ofile)
      if( iprint.ge.1 ) then
        write(*,*) 'Opened output file>', ofile(1:lnblnk(ofile))
      endif
      open(unit=1,file=ofile,status='unknown')

c...write comment block
      write(1,'(a)') '# edi2j'
      write(1,'(a)') '# date: '//date_string
      write(1,'(a)') '#'
      write(1,'(a)') '# file: '//ifile(1:lnblnk(ifile))
      write(1,'(a)') '#'
      write(1,'(a,f6.1)') '# Data in angle of         :', datang
      write(1,'(a,f6.1)') '# Data required in angle of:', angle
      write(1,'(a,f6.1)') '# Rotation applied         :', delta
      write(1,'(a)') '#'
      write(1,'(a,e12.6)') '# MISDAT                  :', misdat
      write(1,'(a)') '#'
      write(1,'(a)') '# INFO block from edi file'
      write(1,'(a)') '#'
      open(unit=77,file=ifile,status="OLD")
      read(77,'(a)') line
      do while( index(line,'>INFO').eq.0 )
        read(77,'(a)') line
      enddo
      read(77,'(a)') line
      do while( index(line(1:1),'>').eq.0 )
        write(1,'(a)') '#'//line(1:lnblnk(line))
        read(77,'(a)') line
      enddo
      write(1,'(a)') '#'

c...write info block
        write(1,'(a  )') '>STATION   ='//site
        write(1,'(a12,f12.1)') '>AZIMUTH   =', angle

c...read latitude after >=DEFINEMEAS line
        rewind(77)
        do while( index(line,'>=DEFINEMEAS').eq.0 )
          read(77,'(a)',end=8) line
        enddo
        goto 9
8       write(*,*)'Cannot find >=DEFINEMEAS line in EDI file'
        stop
9       if( iprint.ge.2 ) then
          write(*,*)'>=DEFINEMEAS line found in EDI file'
        endif
        do while( index(line,'REFLAT=').eq.0 )
          read(77,'(a)',end=10) line
        enddo
        goto 12
10      write(*,*)'REFLAT not found, looking for LAT'
        rewind(77)
        do while( index(line,'LAT').eq.0 )
          read(77,'(a)',end=11) line
        enddo
        goto 12
11      write(*,*)'LAT not found. Setting to 00:00:00'
        line = '00:00:00'
        goto 13
12      line = line(index(line,'=')+1:)
c...drop leading blanks
        do while( line(1:1).eq.' ' )
          line = line(2:)
        enddo
        if(iprint.ge.2) then
          write(*,*)'Latitude line=>',line(1:lnblnk(line))
        endif
13      continue
        xlat = ddmmss2deg(line)
        if(iprint.ge.1) write(*,*)'xlat       =>', xlat
c...read longitude AFTER >=DEFINEMEAS line
        rewind(77)
        do while( index(line,'>=DEFINEMEAS').eq.0 )
          read(77,'(a)',end=8) line
        enddo
        do while( index(line,'REFLONG').eq.0 )
          read(77,'(a)',end=14) line
        enddo
        goto 16
14      write(*,*)'REFLON not found, looking for LON'
        rewind(77)
        do while( index(line,'LON').eq.0 )
          read(77,'(a)',end=15) line
        enddo
        goto 16
15      write(*,*)'LON not found. Setting to 00:00:00'
        line = '00:00:00'
        goto 17
16      line = line(index(line,'=')+1:)
        if(iprint.ge.2) then
          write(*,*)'Longitude line=>',line(1:lnblnk(line))
        endif
17      continue
        xlon = ddmmss2deg(line)
        if(iprint.ge.1) write(*,*)'xlon =>', xlon
c...read elevation
        do while( index(line,'REFELEV').eq.0 )
          read(77,'(a)',end=20) line
        enddo
        leq = index(line,'=')
        read(line(leq+1:),*) elev
        if(iprint.ge.1) write(*,*)'elev =>', elev
        goto 21
20      write(*,*)'REFELEV not found. Set to 0.0'
        elev = 0.0
21      continue

        write(1,'(a12,f12.4)') '>LATITUDE  =', xlat
        write(1,'(a12,f12.4)') '>LONGITUDE =', xlon
        write(1,'(a12,f12.1)') '>ELEVATION =', elev


        write(1,'(a6,1x,f5.1)') site, angle

c-- loop over frequencies 
        itz = 0
        nrhomisdat = 0
        ntzmisdat = 0
        do 100 i = 1, nfreq

          if( SPECTRA ) then
            do j = 1, 7
              do k = 1, 7
                temp(j,k) = sdata(j,k,i)
              enddo
            enddo

            avgs =  info(i,3)
            call mtcomp(temp,freq(i),delta,avgs,Z,vimp,rho,phs,
     &                  srh,sph)
 
c...convert impedances and variances from field units to S.I.
            Z(1,1)    = factor * Z(1,1)
            Z(1,2)    = factor * Z(1,2)
            Z(2,1)    = factor * Z(2,1)
            Z(2,2)    = factor * Z(2,2)
            vimp(1,1) = factor * factor * vimp(1,1)
            vimp(1,2) = factor * factor * vimp(1,2)
            vimp(2,1) = factor * factor * vimp(2,1)
            vimp(2,2) = factor * factor * vimp(2,2)

          else
c...check for misdat
             if( real(zdata(1,i)).eq.misdat .or. 
     &           imag(zdata(1,i)).eq.misdat .or. 
     &           real(zdata(2,i)).eq.misdat .or. 
     &           imag(zdata(2,i)).eq.misdat .or. 
     &           real(zdata(3,i)).eq.misdat .or. 
     &           imag(zdata(3,i)).eq.misdat .or. 
     &           real(zdata(4,i)).eq.misdat .or. 
     &           imag(zdata(4,i)).eq.misdat  ) then
              imiss(i) = 0
              nrhomisdat = nrhomisdat+1
              write(*,*)'MISDAT at freq ', freq(i)
              goto 101
            endif
            imiss(i) = 1

c...convert impedances from field units to S.I.
            if( factor.ne.1. ) then
              zdata(1,i) = factor * zdata(1,i)
              zdata(2,i) = factor * zdata(2,i)
              zdata(3,i) = factor * zdata(3,i)
              zdata(4,i) = factor * zdata(4,i)
              zerr (1,i) = factor * sqrt(zvar(1,i))
              zerr (2,i) = factor * sqrt(zvar(2,i))
              zerr (3,i) = factor * sqrt(zvar(3,i))
              zerr (4,i) = factor * sqrt(zvar(4,i))
            else
              zerr (1,i) = sqrt(zvar(1,i))
              zerr (2,i) = sqrt(zvar(2,i))
              zerr (3,i) = sqrt(zvar(3,i))
              zerr (4,i) = sqrt(zvar(4,i))
            endif

c...write MTs into z array
            if( EGBERT ) then
              Z(1,1) = -zdata(1,i)
              Z(1,2) = -zdata(2,i)
              Z(2,1) = -zdata(3,i)
              Z(2,2) = -zdata(4,i)
            else
              Z(1,1) =  zdata(1,i)
              Z(1,2) =  zdata(2,i)
              Z(2,1) =  zdata(3,i)
              Z(2,2) =  zdata(4,i)
            endif  

c...write GTFs into z array
            Z(3,1) =    zdata(5,i)
            Z(3,2) =    zdata(6,i)

            rho(1,1) = convz2r(Z(1,1), 1./freq(i))
            rho(1,2) = convz2r(Z(1,2), 1./freq(i))
            rho(2,1) = convz2r(Z(2,1), 1./freq(i))
            rho(2,2) = convz2r(Z(2,2), 1./freq(i))
            phs(1,1) = convz2p(Z(1,1), 1./freq(i))
            phs(1,2) = convz2p(Z(1,2), 1./freq(i))
            phs(2,1) = convz2p(Z(2,1), 1./freq(i))
            phs(2,2) = convz2p(Z(2,2), 1./freq(i))
          endif

          rxx(i) = rho(1,1)
          rxy(i) = rho(1,2)
          ryx(i) = rho(2,1)
          ryy(i) = rho(2,2)

          pxx(i) = phs(1,1)
          pxy(i) = phs(1,2)
          pyx(i) = phs(2,1)
          pyy(i) = phs(2,2)

          zxxr(i) = real(Z(1,1))
          zxyr(i) = real(Z(1,2))
          zyxr(i) = real(Z(2,1))
          zyyr(i) = real(Z(2,2))

          zxxi(i) = imag(Z(1,1))
          zxyi(i) = imag(Z(1,2))
          zyxi(i) = imag(Z(2,1))
          zyyi(i) = imag(Z(2,2))

          if( SPECTRA ) then
            if (srh(1,1) .gt. 15 ) srh(1,1) = 15.0
            srxx(i) = srh(1,1)
            if (srh(1,2) .gt. 15 ) srh(1,2) = 15.0
            srxy(i) = srh(1,2)
            if (srh(2,1) .gt. 15 ) srh(2,1) = 15.0
            sryx(i) = srh(2,1)
            if (srh(2,2) .gt. 15 ) srh(2,2) = 15.0
            sryy(i) = srh(2,2)

            if (sph(1,1) .gt. 180 ) sph(1,1) = 180.0
            spxx(i) = sph(1,1)
            if (sph(1,2) .gt. 180 ) sph(1,2) = 180.0
            spxy(i) = sph(1,2)
            if (sph(2,1) .gt. 180 ) sph(2,1) = 180.0
            spyx(i) = sph(2,1)
            if (sph(2,2) .gt. 180 ) sph(2,2) = 180.0
            spyy(i) = sph(2,2)

            vzxx(i) = vimp(1,1)
            vzxy(i) = vimp(1,2)
            vzyx(i) = vimp(2,1)
            vzyy(i) = vimp(2,2)
          else
            vzxx(i) = factor * factor * zvar(1,i)
            vzxy(i) = factor * factor * zvar(2,i)
            vzyx(i) = factor * factor * zvar(3,i)
            vzyy(i) = factor * factor * zvar(4,i)
          endif

c...check for non-zero values of TZ
101       continue
          if( real(Z(3,1)).ne.0. .or. imag(Z(3,1)).ne.0. .or.
     &        real(Z(3,2)).ne.0. .or. imag(Z(3,2)).ne.0. ) then
            itz = itz + 1
            tzxr(itz) = real(Z(3,1))
            tzxi(itz) = imag(Z(3,1))
            tzyr(itz) = real(Z(3,2))
            tzyi(itz) = imag(Z(3,2))
            if( SPECTRA ) then
              vtzx(itz) = vimp(3,1)
              vtzy(itz) = vimp(3,2)
            else
              vtzx(itz) = zvar(5,i)
              vtzy(itz) = zvar(6,i)
            endif
          endif
  
100     continue

c-------------------------------------------------------------------------------
c--write Jones format file

c...change freqs to periods, setting freqs>0 as negative
      do i = 1, nfreq
        if (freq(i) .lt. 1.0 ) then
          freq(i) = 1.0/freq(i)
        else
          freq(i) = -freq(i)
        endif
      enddo

c...arrange in increasing period order
      if( freq(1).lt.freq(nfreq) ) then
        istart = 1
        iend   = nfreq
        iinc   = 1
      else
        istart = nfreq
        iend   = 1
        iinc   = -1
      endif

c...RXX
        write(1,'(A3)') 'RXX'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do 201 i = istart, iend, iinc
           if( SPECTRA ) then
             write(1,1000) freq(i),rxx(i),pxx(i),
     &       10.0**(alog10(rxx(i))+srxx(i)),
     &       10.0**(alog10(rxx(i))-srxx(i)),
     &       pxx(i)+spxx(i), pxx(i)-spxx(i), one, one
           else
             if( imiss(i).eq.0  ) goto 201
             zz    = zdata(1,i)
             zzerr = zerr(1,i)
             if(freq(i).lt.0.) then
               period = -1./freq(i)
             else
               period=freq(i)
             endif
             rhoa  = rxx(i)
             phase = pxx(i)
             zr    = zzerr*cosd(phase)
             zi    = zzerr*sind(phase)
             rhomax = convz2r( zz+cmplx(zr,zi), period )
             rhomin = convz2r( zz-cmplx(zr,zi), period )
             phaerr = atand( zzerr/cabs(zz) )
             phamax = phase + phaerr
             phamin = phase - phaerr
             write(1,1000) freq(i),rhoa,phase, 
     &       rhomax, rhomin, phamax, phamin, one, one
           endif
201      continue

c...RXY
        write(1,'(A3)') 'RXY'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do 202 i = istart, iend, iinc
           if( SPECTRA ) then
             write(1,1000) freq(i),rxy(i),pxy(i), 
     &       10.0**(alog10(rxy(i))+srxy(i)),
     &       10.0**(alog10(rxy(i))-srxy(i)),
     &       pxy(i)+spxy(i), pxy(i)-spxy(i), one, one
           else
             if( imiss(i).eq.0  ) goto 202
             zz    = zdata(2,i)
             zzerr = zerr(2,i)
             if(freq(i).lt.0.) then
               period = -1./freq(i)
             else
               period=freq(i)
             endif
             rhoa  = rxy(i)
             phase = pxy(i)
             zr    = zzerr*cosd(phase)
             zi    = zzerr*sind(phase)
             rhomax = convz2r( zz+cmplx(zr,zi), period )
             rhomin = convz2r( zz-cmplx(zr,zi), period )
             phaerr = atand( zzerr/cabs(zz) )
             phamax = phase + phaerr
             phamin = phase - phaerr
             write(1,1000) freq(i),rhoa,phase, 
     &       rhomax, rhomin, phamax, phamin, one, one
           endif
202      continue

c...RYX
        write(1,'(A3)') 'RYX'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do 203 i = istart, iend, iinc
           if( SPECTRA ) then
             write(1,1000) freq(i),ryx(i),pyx(i),
     &       10.0**(alog10(ryx(i))+sryx(i)),
     &       10.0**(alog10(ryx(i))-sryx(i)),
     &       pyx(i)+spyx(i), pyx(i)-spyx(i), one, one
           else
             if( imiss(i).eq.0  ) goto 203
             zz    = zdata(3,i)
             zzerr = zerr(3,i)
             if(freq(i).lt.0.) then
               period = -1./freq(i)
             else
               period=freq(i)
             endif
             rhoa  = ryx(i)
             phase = pyx(i)
             zr    = zzerr*cosd(phase)
             zi    = zzerr*sind(phase)
             rhomax = convz2r( zz+cmplx(zr,zi), period )
             rhomin = convz2r( zz-cmplx(zr,zi), period )
             phaerr = atand( zzerr/cabs(zz) )
             phamax = phase + phaerr
             phamin = phase - phaerr
             write(1,1000) freq(i),rhoa,phase, 
     &       rhomax, rhomin, phamax, phamin, one, one
           endif
203      continue

c...RYY
        write(1,'(A3)') 'RYY'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do 204 i = istart, iend, iinc
           if( SPECTRA ) then
             write(1,1000) freq(i),ryy(i),pyy(i),
     &       10.0**(alog10(ryy(i))+sryy(i)),
     &       10.0**(alog10(ryy(i))-sryy(i)),
     &       pyy(i)+spyy(i), pyy(i)-spyy(i), one, one
           else
             if( imiss(i).eq.0  ) goto 204
             zz    = zdata(4,i)
             zzerr = zerr(4,i)
             if(freq(i).lt.0.) then
               period = -1./freq(i)
             else
               period=freq(i)
             endif
             rhoa  = ryy(i)
             phase = pyy(i)
             zr    = zzerr*cosd(phase)
             zi    = zzerr*sind(phase)
             rhomax = convz2r( zz+cmplx(zr,zi), period )
             rhomin = convz2r( zz-cmplx(zr,zi), period )
             phaerr = atand( zzerr/cabs(zz) )
             phamax = phase + phaerr
             phamin = phase - phaerr
             write(1,1000) freq(i),rhoa,phase, 
     &       rhomax, rhomin, phamax, phamin, one, one
           endif
204      continue

c...ZXX
        write(1,'(A3)') 'ZXX'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do i = istart, iend, iinc
          if( imiss(i).ne.0  ) then
            write(1,1001) freq(i),zxxr(i),zxxi(i),sqrt(vzxx(i)),one
          endif
        enddo

c...ZXY
        write(1,'(A3)') 'ZXY'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do i = istart, iend, iinc
          if( imiss(i).ne.0  ) then
             write(1,1001) freq(i),zxyr(i),zxyi(i),sqrt(vzxy(i)),one 
          endif
        enddo

c...ZYX
        write(1,'(A3)') 'ZYX'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do i = istart, iend, iinc
          if( imiss(i).ne.0  ) then
             write(1,1001) freq(i),zyxr(i),zyxi(i),sqrt(vzyx(i)),one 
          endif
        enddo

c...ZYY
        write(1,'(A3)') 'ZYY'
        write(1,'(i3.3)') nfreq - nrhomisdat

        do i = istart, iend, iinc
          if( imiss(i).ne.0  ) then
             write(1,1001) freq(i),zyyr(i),zyyi(i),sqrt(vzyy(i)),one 
          endif
        enddo

c...THIS NEEDS MORE WORK IN CASE OF TZ AT SOME FREQS BUT NOT ALL!!!
        if( itz.gt.0 ) then
          if( itz.ne.nfreq) then
            write(*,*)'See Alan - not going to write out TZs properly!'
          endif
c...TZX
          write(1,'(A3)') 'TZX'
          write(1,'(i3.3)') itz 
          do i = istart, iend, iinc
             write(1,1001) freq(i),tzxr(i),tzxi(i),sqrt(vtzx(i)),one 
          enddo

c...TZY
          write(1,'(A3)') 'TZY'
          write(1,'(i3.3)') nfreq 
          do i = istart, iend, iinc
             write(1,1001) freq(i),tzyr(i),tzyi(i),sqrt(vtzy(i)),one 
          enddo
  
        else
          write(*,*)'No TZ estimates in EDI file'
        endif

        if( ONEEDI ) goto 999

1000    format(g11.5,1x,g11.5,1x,f7.2,1x,g11.5,1x,
     &         g11.5,1x,f7.2,1x,f7.2,1x,f4.2,1x,f4.2)
1001    format(g11.5,1x,g11.5,1x,g11.5,1x,
     &         g11.5,1x,f4.2)

        enddo

999     stop
        end




           

