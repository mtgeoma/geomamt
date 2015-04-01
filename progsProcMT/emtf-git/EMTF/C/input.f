      subroutine ininit(inunit)

c        initializes input file: reads header, passes header data
c        to ouput initialization, computes clock zero etc;
c       can be modified for various input formats

      character*7 stname
      character*80 cfin,cfout,chead,cfplt

c  common block HEADER contains info passed to output file, but not
c    used in main program
      character*80 comment
      integer iyr,imo,iday,ihr,imin,isec,iclk0
      real dr
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment

      integer mday(12)
      data mday/31,28,31,30,31,30,31,31,30,31,30,31/
      data stname /'       '/

c**** common block FORMBLK contains info passed to input routine,
c   describing format of input file
      integer itrec
      character*80 cform,cformat
      common /FORMBLK/cform,cformat,itrec,nskip

c  set format 
c      cformat='(5x,3i5)'
      cformat='(5i7)'
c    nskip is number of header records to skip in input data file
      nskip = 0

c*******************************************************************

      print*,'enter input file name'
      read(5,'(a80)') cfin
      stname = cfin(1:6)
      cform = 'formatted'
      open(unit = inunit,file = cfin,form=cform,status = 'old')

c>>>>>>>>>>>>>>>>    read header info from h_ file
      print*,'enter name for clock reset file'
      read(5,'(a80)') chead
      open(unit = 10, file = chead)
c*nsampling rate (seconds)
      read(10,*) dr
c*instrument clock zero time - year month, day, time (ut) [ hour,min,sec]
      read(10,*) iyr,imo,iday,ihr,imin
c*universal clock zero time - year month, day, time (ut) [ hour,min,sec]
      read(10,*) iyru,imou,idayu,ihru,iminu
      isec = 0
      isecu = 0
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c     compute julian day number for instrument clock zero
        jday = 0
        do 20 i = 1,imo-1
        if( (mod(iyr,4).eq.0) .and. (i.eq.2) ) then
           jday = jday+mday(i) + 1
        else
           jday = jday + mday(i)
        end if
20      continue
        jday = jday + iday

c     compute julian day number for universal clock zero
        jdayu = 0
        do 25 i = 1,imou-1
        if( (mod(iyru,4).eq.0) .and. (i.eq.2) ) then
           jdayu = jdayu+mday(i) + 1
        else
           jdayu = jdayu + mday(i)
        end if
25      continue
        jdayu = jdayu + idayu

c   record number (relative to universal clock zero) of instrument
c      clock zero

      jday = jday-jdayu + 1
      ihr = ihr - ihru
      imin = imin - iminu
      isec = isec - isecu
      iclk0 = nint((((((jday-1)*24+ihr)*60+imin)*60+isec))/dr)
      close(10)

c           skip header records
      do 30 i = 1,nskip
30          read(inunit,*)
c       use itrec for keeping track of record # for files with
c        no record #s
      itrec = iclk0-1
      return
      end
c_____________________________________________________________________c
      subroutine frstdat(inunit)
   
c        position input file at firstheader record

      character*80 comment
      integer iyr,imo,iday,ihr,imin,isec,iclk0
      real dr
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment

      integer itrec
      character*80 cform,cformat
      common /FORMBLK/cform,cformat,itrec,nskip

      rewind(inunit)

      do 30 i = 1,nskip
30       read(inunit,*)

c       use itrec for keeping track of record # for files with
c        no record #s
      itrec = iclk0-1
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

      integer itrec
      character*80 cform,cformat
      common /FORMBLK/cform,cformat,itrec,nskip
        
      nch=nchp1-1
      lend = .false.
      j = ipoint
1     continue
      read(inunit,cformat,end=100) (ix1(k,j),k=2,nchp1)
      itrec = itrec + 1
      ix1(1,j) = itrec
      j = j+1
      if(mod(itrec,120).eq.0) then
         read(inunit,cformat) idum1,idum2,idum3
      endif
      if(j.le.n) go to 1
      ngot = j - 1
ccc       n will be 120 (so 10 minuites worth of data are read
ccc       skip the next record (time mark)
ccc       (for now print out to see if things work)
      return
        
100   ngot = j - 1 
      lend = .true.
      return
      end
