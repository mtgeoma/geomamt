      subroutine ediresread(ind,freq,nfreq,zdata,zvar,info,site)
c
c    Subroutine to read a response EDI file
c    Alan Jones, 23 January 2001
c---------------------------------------------------------------------

      implicit none

c---------------------------------------------------------------------

      include 'ctrlblk.inc'
      include 'size.inc'

c---------------------------------------------------------------------

      character line*80, site*6
      character dataid*10, elev*10, lat*15, lon*15
      real freq(MAXDAT), info(MAXDAT,3)
      real datar(MAXDAT),datai(MAXDAT)
      complex zdata(6,MAXDAT) 
      real zvar(6,MAXDAT)
      integer i, nfreq, ind

c     external lnblnk
      integer lnblnk

      logical search/.TRUE./
      logical DEBUG

      common /cmndbg/ DEBUG

c==========================================================
c-- find and read header

      if( DEBUG ) write(*,*)'ediresread: entered'
       
      rewind( ind )

c-- find and read NFREQ

      search = .TRUE.
      do while(search)

        read(ind,'(a80)') line
        call shiftup(line)

        if (index(line,'NFREQ=').gt.0) then
c...read NFREQ
          read(line(7:),*) nfreq
          search = .FALSE.
        endif

      enddo
      if(iprint.ge.1) write(*,*)'nfreq =', nfreq
      backspace(ind)

c-- find and read FREQ

      search = .TRUE.
      do while(search)

        read(ind,'(a80)') line
        call shiftup(line)
        if( iprint.ge.3 ) then
          write(*,*)'line=>', line(1:lnblnk(line))
          write(*,*)'lnblnk(line)=>',lnblnk(line)
        endif

        if (index(line,'>ZROT').gt.0) search=.FALSE.

        if (index(line,'>FREQ').gt.0) then
c...read NFREQ
          read(ind,*) (freq(i),i=1,nfreq)
	  if( iprint.ge.1 ) then
            write(*,*)'Read in freqs ==>'
            write(*,*) (freq(i),i=1,nfreq)
	  endif  
          search = .FALSE.
        endif

      enddo
      backspace(ind)

c-- find and read the transfer functions

      search = .TRUE.
      do while(search)

        read(ind,'(a80)',end=999) line
        call shiftup(line)
        if( iprint.ge.3 ) then
          write(*,*)'line=>', line(1:lnblnk(line))
          write(*,*)'lnblnk(line)=>',lnblnk(line)
        endif

        if (index(line,'>END').gt.0) search=.FALSE.

c...ZXX
        if (index(line,'>ZXXR').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZXXR...'
          read(ind,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZXXI').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZXXI...'
          read(ind,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(1,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZXX.VAR').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZXX.VAR...'
          read(ind,*) (zvar(1,i),i=1,nfreq)
        endif

c...ZXY
        if (index(line,'>ZXYR').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZXYR...'
          read(ind,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZXYI').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in ZXYI...'
          read(ind,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(2,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZXY.VAR').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZXY.VAR...'
          read(ind,*) (zvar(2,i),i=1,nfreq)
        endif

c...ZYX
        if (index(line,'>ZYXR').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZYXR...'
          read(ind,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZYXI').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in ZYXI...'
          read(ind,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(3,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZYX.VAR').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZYX.VAR...'
          read(ind,*) (zvar(3,i),i=1,nfreq)
        endif

c...ZYY
        if (index(line,'>ZYYR').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in ZYYR...'
          read(ind,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>ZYYI').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in ZYYI...'
          read(ind,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(4,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>ZYY.VAR').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZYY.VAR...'
          read(ind,*) (zvar(4,i),i=1,nfreq)
        endif

c...TZX
        if (index(line,'>TXR.EXP').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in TXR...'
          read(ind,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>TXI.EXP').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in TXI...'
          read(ind,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(5,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>TXRVAR.EXP').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZXRVAR...'
          read(ind,*) (zvar(5,i),i=1,nfreq)
        endif

c...TZY
        if (index(line,'>TYR.EXP').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in TYR...'
          read(ind,*) (datar(i),i=1,nfreq)
        endif
        if (index(line,'>TYI.EXP').gt.0) then
          if( iprint.ge.1 )write(*,*)'Reading in TYI...'
          read(ind,*) (datai(i),i=1,nfreq)
          do i = 1, nfreq
            zdata(6,i) = cmplx( datar(i), datai(i) )
          enddo
        endif
        if (index(line,'>TYRVAR.EXP').gt.0) then
           if( iprint.ge.1 )write(*,*)'Reading in ZYRVAR...'
          read(ind,*) (zvar(6,i),i=1,nfreq)
        endif

      enddo

999   close(ind)   
      if( DEBUG ) write(*,*)'ediresread: exiting'
      return
      end
