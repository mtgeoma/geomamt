c**********************************************************
c
      subroutine edispecread(ifile,misdat,freq,nfreq,sdata,info,site)
c
c    Subroutine to read a Spectral EDI file
c    Gary McNeice July, 1996
c**********************************************************

c... if you wish to change MAXF here, you also have to change it in the main code
      integer MAXF
      parameter ( MAXF = 200 )

      character line*80, ifile*80, site*6, temp*20
      character dataid*20, elev*10, lat*10, lon*10
      character maxchan*10, char1*1
      real freq(MAXF), info(MAXF,3)
      real sdata(7,7,MAXF) 
      real misdat
      integer nlen, nchan, count, j, nfreq, ltemp
      logical search/.true./
      
      common /ctrl_blk/ iprint

c==========================================================

c-- open EDI file

      open(unit=77,file=ifile,status="OLD")
      if( iprint.ge.1 ) write(*,*)'edispecread: opened edi file>',
     &                             ifile(1:lnblnk(ifile))

c-- find and read header

      do while(search)

        read(77,'(a80)') line
        call shiftup(line)

        if( index(line,'EMPTY').gt.0 ) then
          read(line(index(line,'=')+1:),*) misdat
        endif

        if (index(line,'>INFO').gt.0) search=.false.

        if (index(line,'DATAID').gt.0) then
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

          if( iprint.ge.3 ) write(*,*) 'edispecread: DataID = ', 
     &                                 dataid(1:lnblnk(dataid))
          site(1:3) = dataid(1:3)
          call shiftdown(site)
          if( iprint.ge.3 ) write(*,*) 'edispecread: site = ', site(1:3)
        endif

        if (index(line,'LAT').gt.0) then
          lat = line(index(line,'=')+1:len(line))
          call trunc(lat,nlen)
        endif

        if (index(line,'LONG').gt.0) then
          lon = line(index(line,'=')+1:len(line))
          call trunc(lon,nlen)
        endif

        if (index(line,'ELEV').gt.0) then
          elev = line(index(line,'=')+1:len(line))
        end if

      enddo
      if( iprint.ge.1 ) then
        write(*,*)'edispecread: lat  = ', lat
        write(*,*)'edispecread: lon  = ', lon
        write(*,*)'edispecread: elev = ', elev
      endif
      backspace(77)

       
c-- find and read MAXCHAN from DEFINEMEAS
      read(77,'(a80)') line
      call shiftup(line)
      do while( index(line,'>=DEFINEMEAS').eq.0 )
        read(77,'(a80)',end=9902) line
        call shiftup(line)
      enddo
      if( iprint.ge.2 ) write(*,*)'edispecread: found DEFINEMEAS'
      
      read(77,'(a80)') line
      call shiftup(line)
      do while( index(line,'MAXCHAN=').eq.0 )
        read(77,'(a80)') line
        call shiftup(line)
      enddo
      if( iprint.ge.2 ) write(*,*)'edispecread: found MAXCHAN: ',
     &                            line(1:lnblnk(line))
            
      maxchan = line(index(line,'=')+1:index(line,'=')+1)
      call atoi(maxchan,nchan)
      if( iprint.ge.1 ) write(*,*)'edispecread: MAXCHAN, NCHAN =',
     &                             maxchan, nchan
      if (nchan.gt.7 ) then
        write(*,*) 'More than 7x7 spectral matrix found !'
        write(*,*) 'ERROR READING MATRIX !'
        stop
      endif
      if( nchan.le.0 ) then
        write(*,*) 'edispecread: Error nchan = 0'
        stop
      endif

      rewind(77)

c-- find and read HMEAS and EMEAS

      search = .true.
      do while(search)

        read(77,'(a80)') line
        call shiftup(line)

        if (index(line,'>=SPEC').gt.0) search=.false.
        
        if (index(line,'>HMEAS').gt.0) then
          site(4:6) = line(index(line,'.')+1:index(line,'.')+3)
          search = .false.
        endif
      
      enddo
      rewind(77)

      write(*,*)'edispecread: site = ', site

c-- find and read NFREQ from SPECTRASECT
c...find start of SPECTRASECT
      read(77,'(a80)') line
      call shiftup(line)
      do while( index(line,'>=SPECTRASECT').eq.0 )
        read(77,'(a80)',end=9901) line
        call shiftup(line)
      enddo
      if( iprint.ge.2 ) write(*,*)'edispecread: found SPECTRASECT'

c...skip over SECTID
      read(77,'(a80)') line
c...skip over NCHAN
      read(77,'(a80)') line
c...read in NFREQ line
      read(77,'(a80)') line
      call shiftup(line)
      if( iprint.ge.5 ) write(*,*)'NFREQ line>', line(1:lnblnk(line))
      temp = line(index(line,'NFREQ=')+6:)
c...drop blanks at beginning
      do while( temp(1:1).eq.' ' )
        temp = temp(2:)
      enddo
      call atoi(temp,nfreq)
      if( iprint.ge.2 ) write(*,*)'edispecread: NFREQ =', nfreq

      backspace(77)

c-- find and read the SPECTRA

      count = 0
      search = .true.
      do while(search)

        read(77,'(a80)',end=999) line
        call shiftup(line)

        if (index(line,'>END').gt.0) search=.false.

        if (index(line,'>SPECTRA').gt.0) then
          count = count + 1
          if( iprint.ge.5 ) write(*,*)'edispecread: count =', count
c...get frequency
          temp = line(index(line,'FREQ=')+5:)
c...get rid of blank spaces at beginning
          do while( temp(1:1).eq.' ') 
            temp = temp(2:)
          enddo
          temp = temp(1:index(temp,' ')-1)
          read(temp(1:lnblnk(temp)),*) freq(count)
          if( iprint.ge.5 ) write(*,*)'edispecread: FREQ =', freq(count)

c...get ROTSPEC
          if( index(line,'ROTSPEC=').eq.0 ) then
            if( count.eq.1 ) then
              write(*,*)'edispecread: ROTSPEC not found - set to zero'
            endif
            info(count,1) = 0
          else
            temp = line(index(line,'ROTSPEC=')+8:)
c...get rid of blank spaces at beginning
            do while( temp(1:1).eq.' ') 
              temp = temp(2:)
            enddo  
            ltemp = index(temp,' ')-1
            read(temp(1:ltemp),*) info(count,1)
c          call atof(temp,info(count,1))
            if( iprint.ge.5 ) write(*,*)'edispecread: ROTSPEC =', 
     &                                info(count,1)
          endif

c...BW is the frequency bandwidth
          temp = line(index(line,'BW=')+3:)
c...get rid of blank spaces at beginning
          do while( temp(1:1).eq.' ') 
            temp = temp(2:)
          enddo
          ltemp = index(temp,' ')-1
          read(temp(1:ltemp),*) info(count,2)
c          call atof(temp,info(count,2))
          if( iprint.ge.3 ) write(*,*)'edispecread: BW =', 
     &                                info(count,2)

c...AVGT is the number of degrees of freedom in the spectral estimate. If this is missing, then the variances cannot be calculated
          temp = line(index(line,'AVGT=')+5:)
c...get rid of blank spaces at beginning
          do while( temp(1:1).eq.' ') 
            temp = temp(2:)
          enddo
          ltemp = index(temp,' ')-1
          read(temp(1:ltemp),*) info(count,3)
c          call atof(temp,info(count,3))
          if( iprint.ge.3 ) write(*,*)'edispecread: AVGT =', 
     &                                info(count,3)

          read(77,*) ((sdata(i,j,count),j=1,nchan),i=1,nchan)
c         do i = 1, nchan 
c            read(77,*) (sdata(i,j,count),j=1,nchan)
c            read(77,666) (sdata(i,j,count),j=1,7)
c666         format(1X,E10.3,1X,E10.3,1X,E10.3,1X,E10.3,1X,E10.3,
c     &             1X,E10.3,1X,E10.3)
c          enddo    
        endif
      enddo

999   nfreq = count
      close(77)     
      return
      
9901  write(*,*)'***Could not find SPECTRASECT in edi file***'
      stop
      
9902  write(*,*)'***Could not find DEFINEMEAS in edi file***'
      stop
      end
     
c--------------------------------------------------------------------------
      subroutine atof(string,value)
c--------------------------------------------------------------------------
      character string*(*), dumm*20, fmt*8
      real      value
      integer   nc, nd
c
      dumm = string
      nc = index(dumm,' ') - 1
      if(nc.le.0) call trunc(dumm,nc)
      dumm = dumm(1:nc)
      nd = index(dumm,'.')
      if(nd.ne.0) nd = nc - nd
      write(fmt,1) nc, nd
1     format('(F',i2,'.',i2,')')
      read(dumm(1:nc),fmt,err=1000) value
      return
1000  write(*,*) 'ATOF - Error reading floating point number from ',
     >           dumm(1:nc)
      return
      end
c--------------------------------------------------------------------------
      subroutine atoi(string,ivalue)
c--------------------------------------------------------------------------
      character string*(*), dumm*20, fmt*5
      integer   nc, ivalue
c
      dumm = string
      nc = index(dumm,' ') - 1
      if(nc.le.0) call trunc(dumm,nc)
      dumm = dumm(1:nc)
      write(fmt,1) nc
1     format('(I',i2,')')
      read(dumm(1:nc),fmt,err=1000) ivalue
      return
1000  write(*,*) 'ATOI - Error reading integer from ',dumm(1:nc)
      return
      end

