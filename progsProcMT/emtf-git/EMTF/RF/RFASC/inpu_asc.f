ccc     THIS VERSION of INPUT :::  reads an ASCII file
ccc     with format = cformat (set in format.h) ... assumes one
ccc     sample (all channels in integers) per line, no record number,
ccc     all data points contiguous (no gaps, no time markers),
ccc     no data file header.

ccc     A separate "clock reset file"  is required.  The clock
ccc     file has three lines: (1) sampling rate in seconds (NOT hz !!!!),
ccc     (2) the time of the first record in the file (yr,mo,day,hr,min, all integers)
ccc     (3) the time of record number zero (yr,mo,day,hr,min, all integers)
ccc           (this will be the same for all simultaneous stations, and serves
ccc            to align data from separate files)
ccc     With this reformatting program you have to make system parameter files
ccc     yourself

      subroutine ininit(inunit,cfout,lsp)

c        initializes input file: reads header, passes header data
c        to ouput initialization, computes clock zero etc;
c       can be modified for various input formats
ccc    in this version cfout and lsp are not referred to 


      character*80 cfile,cfhd,cfout
      character*40 cformat

c  common block HEADER contains info passed to output file, but not
c    used in main program
      character*80 comment
      integer iyr,imo,iday,ihr,imin,isec,iclk0

      real dr,clock(2)
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment
      integer mday(12),irec
      data mday/31,28,31,30,31,30,31,31,30,31,30,31/
      common /FORMBLK/irec,cformat

      logical lsp

      include 'format.inc'

1     continue
      print*,'input file name'
      read(5,'(a80)') cfile
      open(unit=inunit,file=cfile,form='formatted',status='old')
 
c>>>>>>>>>>>>>>>>    read header info from clock reset file
      print*,'enter clock reset file name'
      read(5,'(a80)') cfhd
      open(unit = 10, file = cfhd)

ccc   sampling rate (seconds)
      read(10,*) dr

ccc   instrument clock zero time - year month, day, time (ut) [ hour,min,sec]
ccc    (i.e., the time of the first record in the file)
      read(10,*) iyr,imo,iday,ihr,imin

cc    universal clock zero time - year month, day, time (ut) [ hour,min,sec]
ccc    ( the time to call record # 0 )
      read(10,*) iyru,imou,idayu,ihru,iminu

      isec = 0
      isecu = 0
c        now see if an offset (number of records) is in start time file:
      read(10,*,end=25) iroff
c        can also enter clock drift parameters as a final entry
      read(10,*,end=25) clock
25    iroff = 0
      clock(1) = 0.
      clock(2) = 0.
      close(10)
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c     compute julian day number for instrument clock zero
      jday = 0
      do i = 1,imo-1
         if( (mod(iyr,4).eq.0) .and. (i.eq.2) ) then
            jday = jday+mday(i) + 1
         else
            jday = jday + mday(i)
         end if
      enddo
      jday = jday + iday

cc     compute julian day number for universal clock zero
      jdayu = 0
      do i = 1,imou-1
         if( (mod(iyru,4).eq.0) .and. (i.eq.2) ) then
            jdayu = jdayu+mday(i) + 1
         else
            jdayu = jdayu + mday(i)
         end if
      enddo
      jdayu = jdayu + idayu

c   record number (relative to universal clock zero) of first record in file
      jday = jday-jdayu + 1
      ihr = ihr - ihru
      imin = imin - iminu
      isec = isec - isecu
ccc   irec is absolute record number of first sample in data file
      irec = (((((jday-1)*24+ihr)*60+imin)*60+isec))/dr
      irec = irec + iroff
      return
      end
c_____________________________________________________________________c
      subroutine frstdat(inunit)
   
c        position input file at first header record

      character*80 comment
      character*40 cformat
      integer iyr,imo,iday,ihr,imin,isec,iclk0
      real dr
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment

      integer irec
      common /FORMBLK/irec,cformat

c        irelrc is the relative data record number (relative to start
c      of file .... ircrel = 1 for first data point in file)

      rewind(inunit)
c         (this routine is basically a dummy in this implementation)

      return
      end
c_____________________________________________________________________c
      subroutine indo(inunit,ix1,n,ngot,ipoint,nchp1,lend)
c      inunit = input unit number
c      ix1 = data array (integer)
c      n = number of points to try returning in ix1
c      ngot = number of points actually returned
c      nch = number of data channels
c      ipoint = pointer to starting index in array ix1 (smaller
c         indices in array already have data)
c      lend = .true. if EOF is reached (logical variable)
c      

      integer ix1(nchp1,n)
      logical lend

      integer irec
      character*40 cformat
      common /FORMBLK/irec,cformat
        
      lend = .false.
1     continue
      do j = ipoint,n
         read(inunit,cformat,end=20) (ix1(k,j),k=2,nchp1)
         ix1(1,j) = irec
         irec = irec + 1
      enddo
      ngot = n
      return

20    continue
      ngot = j-1
      lend = .true.
      return
      end
