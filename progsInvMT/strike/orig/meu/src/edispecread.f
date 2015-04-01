      subroutine edispecread(ind,freq,nfreq,sdata,info,site)
c
c    Subroutine to read a Sectral EDI file
c    Gary McNeice July, 1996
c**********************************************************

      implicit none

c----------------------------------------------------------

      include 'size.inc'
      include 'ctrlblk.inc'

c----------------------------------------------------------

      character line*80, site*6, temp*10
      character dataid*10, elev*10, lat*10, lon*10
      character maxchan*10
      real freq(MAXDAT), info(MAXDAT,3)
      real sdata(7,7,maxf) 
      integer i, nlen, nchan, count, j, nfreq, ind

      logical search/.TRUE./
      logical DEBUG

      common /cmndbg/ DEBUG


c==========================================================
c-- find and read header

      if( DEBUG ) write(*,*)'edispecread: entered'

      do while(search)

        read(ind,'(a80)') line
        call shiftup(line)

        if (index(line,'>INFO').gt.0) search=.false.

        if (index(line,'DATAID').gt.0) then
          dataid = line(index(line,'=')+1:len(line))
          site(1:3) = dataid(2:4)
        end if

        if (index(line,'LAT').gt.0) then
          lat = line(index(line,'=')+1:len(line))
          call trunc(lat,nlen)
        end if

        if (index(line,'LONG').gt.0) then
          lon = line(index(line,'=')+1:len(line))
          call trunc(lon,nlen)
        end if

        if (index(line,'ELEV').gt.0) then
          elev = line(index(line,'=')+1:len(line))
        end if

      end do
      backspace(ind)

       
c-- find and read definemeas

      search = .TRUE.

      do while(search)

        read(ind,'(a80)') line
        call shiftup(line)

        if (index(line,'>HMEAS').gt.0) search=.false.
        if (index(line,'>EMEAS').gt.0) search=.false.

        if (index(line,'MAXCHAN').gt.0) then
          maxchan = line(index(line,'=')+1:len(line))
          call atoi(maxchan,nchan)
          if (nchan .gt. 7 ) then
             write(*,*) 'More than 7x7 spectral matrix found !'
             write(*,*) 'ERROR READING MATRIX !'
             stop
          end if
        end if

      end do
      backspace(ind)

c-- find and read HMEAS and EMEAS

      search = .true.
      do while(search)

        read(ind,'(a80)') line
        call shiftup(line)

        if (index(line,'>=SPEC').gt.0) search=.false.
        
        if (index(line,'>HMEAS').gt.0) then
          site(4:6) = line(index(line,'ID=')+8:index(line,'ID=')+10)
          search = .false.
        end if
      
      end do
      backspace(ind)


c-- find and read SPECTRASECT

      search = .true.
      do while(search)

        read(ind,'(a80)') line
        call shiftup(line)

        if (index(line,'>SPECTRA').gt.0) search=.false.

        if (index(line,'>=SPEC').gt.0) then
          read(ind,'(a80)') line
          read(ind,'(a80)') line
          read(ind,'(a80)') line
          call shiftup(line)
          temp = line(index(line,'EQ=')+3:index(line,'EQ=')+4)
          call atoi(temp,nfreq)
        end if

      end do
      backspace(ind)

c-- find and read the SPECTRA

      count = 0
      search = .true.
      do while(search)

        read(ind,'(a80)',end=999) line
        call shiftup(line)

        if (index(line,'>END').gt.0) search=.false.

        if (index(line,'>SPEC').gt.0) then
          count = count + 1
          temp = line(index(line,'EQ=')+3:index(line,'ROT')-2)
          call atof(temp,freq(count))
          temp = line(index(line,'EC=')+3:index(line,'BW=')-2)
          call atof(temp,info(count,1))
          temp = line(index(line,'BW=')+3:index(line,'AVG')-2)
          call atof(temp,info(count,2))
          temp = line(index(line,'GT=')+3:index(line,'//')-2)
          call atof(temp,info(count,3))
          
          do i = 1, nchan 
            read(ind,666) (sdata(i,j,count),j=1,7)
666         format(1X,E10.3,1X,E10.3,1X,E10.3,1X,E10.3,1X,E10.3,
     &             1X,E10.3,1X,E10.3)
          end do    
        end if
      end do

999   nfreq = count
      close(ind) 
      if( DEBUG ) write(*,*)'edispecread: exiting'    
      return
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

