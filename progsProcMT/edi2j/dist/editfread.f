c**********************************************************
c
      subroutine editfread(ifile,misdat,freq,nfreq,zdata,zvar,info,site)
c
c    Subroutine to read a response EDI file
c    Alan Jones, 23 January 2001
c
c---------------------------------------------------------------------

      implicit none

c---------------------------------------------------------------------

      include 'ctrlblk.inc'

c---------------------------------------------------------------------

c... if you wish to change MAXF here, you also have to change it in the main code
      integer MAXF
      parameter ( MAXF = 200 )

      character line*80, ifile*80, site*6
      character dataid*20, elev*10, lat*15, lon*15
      character char1*1
      real freq(MAXF), info(MAXF,3)
      real datar(MAXF),datai(MAXF)
      complex zdata(6,MAXF) 
      real zvar(6,MAXF)
      real misdat
      integer i, nfreq

c     external lnblnk
      integer lnblnk

      logical search/.TRUE./
c     logical FINDSITE/.TRUE./

c==========================================================

c-- open EDI file

      open(unit=77,file=ifile,status="OLD")
      if(iprint.ge.1) then
        write(*,*)'Entered editfread and opened file >', ifile
      endif
      dataid = "XXX000"
      site   = "XXX000"

c-- find and read header

      search = .TRUE.
      do while(search)

        read(77,'(a80)') line
        if(iprint.ge.20) write(*,*)'Line>',line(1:lnblnk(line))
        line = line(1:lnblnk(line))
        call shiftup(line)
        if( iprint.ge.20) then
          write(*,*)'lnblnk(line),line=>',lnblnk(line),
     &              line(1:lnblnk(line))
        endif

        if( index(line,'EMPTY').gt.0 ) then
          read(line(index(line,'=')+1:),*) misdat
        endif

        if( index(line,'>=MTSECT').gt.0 ) search=.FALSE.

c        if( index(line,'DATAID').gt.0 .and. FINDSITE ) then
        if( index(line,'DATAID').gt.0  ) then
          dataid = line(index(line,'=')+1:len(line))
c...rid DataID of leading non-alphanumeric characters
          char1 = dataid(1:1)
          do while( ichar(char1).le.47 .or.
     &       (ichar(char1).ge.58 .and. ichar(char1).le.64) .or.
     &       (ichar(char1).ge.91 .and. ichar(char1).le.96) .or.
     &        ichar(char1).ge.123 ) 
            dataid = dataid(2:)
            char1 = dataid(1:1)
          enddo
c...rid DataID of trailing non-alphanumeric characters
          char1 = dataid(lnblnk(dataid):lnblnk(dataid))
          do while( ichar(char1).le.47 .or.
     &       (ichar(char1).ge.58 .and. ichar(char1).le.64) .or.
     &       (ichar(char1).ge.91 .and. ichar(char1).le.96) .or.
     &        ichar(char1).ge.123 ) 
            dataid(lnblnk(dataid):lnblnk(dataid)) = ' '
            char1 = dataid(lnblnk(dataid):lnblnk(dataid))
          enddo

          if( iprint.ge.3 ) write(*,*) 'editfread: DataID = ', 
     &                                 dataid(1:lnblnk(dataid))
          site(1:3) = dataid(1:3)
          site(4:6) = dataid(lnblnk(dataid)-2:)
          call shiftdown(site)
          if( iprint.ge.3 ) write(*,*) 'editfread: site = ', site(1:3)
c          call tin('Give site name (6-char)', site )
c          FINDSITE = .FALSE.
        endif

        if (index(line,'LAT=').gt.0) then
          lat = line(index(line,'=')+1:lnblnk(line))
          if( iprint.ge.1 ) write(*,'(2a)')'Lat  =>', lat
        endif

        if (index(line,'LONG=').gt.0) then
          lon = line(index(line,'=')+1:lnblnk(line))
          if( iprint.ge.1 ) write(*,'(2a)')'Long =>', lon
        endif

        if (index(line,'ELEV=').gt.0) then
          elev = line(index(line,'=')+1:lnblnk(line))
          if( iprint.ge.1 ) write(*,'(2a)')'Elev =>', elev
        endif

      enddo
      if( iprint.ge.1 ) then
        write(*,'(4a)')'DataID = ', dataid, '    Site = ', site
        write(*,'(6a)')'Lat =', lat,'   Lon =', lon, '  Elev =', elev
      endif
      rewind(77)

       

c-- find and read HMEAS and EMEAS

      search = .TRUE.
      do while(search)

        read(77,'(a80)') line
        call shiftup(line)

        if (index(line,'>=MTSECT').gt.0) search=.FALSE.
        
        if (index(line,'>HMEAS').gt.0) then
c          site(4:6) = line(index(line,'.')+1:index(line,'.')+3)
          search = .FALSE.
        endif
      
      enddo
c      write(*,*)'site = ', site
      backspace(77)

      write(*,*)'editfread: site = ', site

c-- find and read NFREQ

      search = .TRUE.
      do while(search)

        read(77,'(a80)') line
        call shiftup(line)

        if (index(line,'NFREQ=').gt.0) then
c...read NFREQ
          read(line(7:),*) nfreq
          search = .FALSE.
        endif
      enddo
      if(iprint.ge.1) write(*,*)'nfreq =', nfreq
      backspace(77)

c-- find and read FREQ

      search = .TRUE.
      do while(search)

        read(77,'(a80)') line
        call shiftup(line)
        if( iprint.ge.20 ) then
          write(*,*)'line=>', line(1:lnblnk(line))
          write(*,*)'lnblnk(line)=>',lnblnk(line)
        endif

        if (index(line,'>FREQ').gt.0) then
c        write(*,*)'Reading in freqs...'
c...read NFREQ
          read(77,*) (freq(i),i=1,nfreq)
          if( iprint.ge.1 ) then
            write(*,*)'Read in freqs ==>'
            write(*,*) (freq(i),i=1,nfreq)
            endif  
          search = .FALSE.
        endif

      enddo


c-- find and read ZROT

      rewind(77)
      search = .TRUE.
      do while(search)

        read(77,'(a80)',end=10) line
        call shiftup(line)
        if( iprint.ge.20 ) then
          write(*,*)'line=>', line(1:lnblnk(line))
          write(*,*)'lnblnk(line)=>',lnblnk(line)
        endif

        if (index(line,'>ZROT').gt.0) then
          read(77,*) (info(i,1),i=1,nfreq)
          if( iprint.ge.1 ) then
            write(*,*)'Read in ZROT ==>'
            write(*,*) (info(i,1),i=1,nfreq)
          endif  
          search = .FALSE.
        endif

      enddo
      goto 20
      
10    write(*,*)'ZROT not found in file. Set to 0.0'
      do i = 1, nfreq
        info(i,1) = 00.
      enddo


c-- find and read the transfer functions

20    rewind(77)
      search = .TRUE.
      do while(search)

        read(77,'(a80)',end=999) line
        call shiftup(line)
        if( iprint.ge.20 ) then
          write(*,*)'line=>', line(1:lnblnk(line))
          write(*,*)'lnblnk(line)=>',lnblnk(line)
        endif

        if (index(line,'>END').gt.0) search=.FALSE.

c...ZXX
        if (index(line,'>ZXXR').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZXXR...'
          read(77,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZXXI').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZXXI...'
          read(77,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(1,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZXX.VAR').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZXX.VAR...'
          read(77,*) (zvar(1,i),i=1,nfreq)
        endif

c...ZXY
        if (index(line,'>ZXYR').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZXYR...'
          read(77,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZXYI').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in ZXYI...'
          read(77,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(2,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZXY.VAR').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZXY.VAR...'
          read(77,*) (zvar(2,i),i=1,nfreq)
        endif

c...ZYX
        if (index(line,'>ZYXR').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZYXR...'
          read(77,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZYXI').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in ZYXI...'
          read(77,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(3,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZYX.VAR').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZYX.VAR...'
          read(77,*) (zvar(3,i),i=1,nfreq)
        endif

c...ZYY
        if (index(line,'>ZYYR').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in ZYYR...'
          read(77,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZYYI').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in ZYYI...'
          read(77,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(4,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZYY.VAR').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in ZYY.VAR...'
          read(77,*) (zvar(4,i),i=1,nfreq)
        endif

c...TZX
        if (index(line,'>TXR.EXP').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in TXR...'
          read(77,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>TXI.EXP').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in TXI...'
          read(77,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(5,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>TXVAR.EXP').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in TXVAR...'
          read(77,*) (zvar(5,i),i=1,nfreq)
        endif

c...TZY
        if (index(line,'>TYR.EXP').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in TYR...'
          read(77,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>TYI.EXP').gt.0) then
          if( iprint.ge.5 )write(*,*)'Reading in TYI...'
          read(77,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(6,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>TYVAR.EXP').gt.0) then
           if( iprint.ge.5 )write(*,*)'Reading in TYVAR...'
          read(77,*) (zvar(6,i),i=1,nfreq)
        endif

      enddo

      if(iprint.ge.1) then
        write(*,*)'editfread: all data read in - returning to main code'
      endif

999   close(77)   
      return
      end
