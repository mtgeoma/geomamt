c______________________________________________________________________
c
      subroutine cininit(msval,rmsval,nch)

      include 'iounits.inc'
      include '../include/datsz.inc'
      include 'input.inc'

      character*80 cfin
      character*40 ctemp,croot,cfpaths
      character*3 cfn
      logical lpath
c      lfmiss is true to fill missing values
c       cfsp, cfbr,cfdecset,cfpwset,cdirout are path/file names for
c       system params, bad recs,decset,pwset, and output file

c       msval (here set equal to maxint) is missing value code
      msval = 2147483647
      rmsval = msval - 10000

      nmsmx = 250
      lfirst = .true.
      ioff = 0
      lfmiss = .false.

      write(6,'(''Enter input file name:  '',$)')
      read(5,'(a40)') croot
c        input file is just local name; rest of path is found from
c       file cf_paths (if cf_paths isn't in local directory, all files
c       are assumed to be in local directory)

c     changed default name for the paths file (more DOS like)
      cfpaths = 'paths.cfg'
      open(unit=pth_unit,file=cfpaths,status='old',err=20)
      lpath = .true.
      goto 21
cnew    prompt the user for an alternate paths file name
 20   print*,'Couldn''t find the file paths.cfg. Please enter a'
      print*,'NEW NAME for the paths file or hit return and the '
      print*,'program writes output into the working directory.' 
      read(5,'(a40)') cfpaths
      if (cfpaths(1:1).eq.' ') then
         lpath = .false.
      else
         open(unit=pth_unit,file=cfpaths,status='old',err=20)
         lpath = .true.
      endif

c        data directory
 21   continue

ccc   mi is length of data file name with blanks stripped off
      mi = iclong(croot,40)
ccc   mr is length of "file root" -- data file name with suffix (beginging
ccc     with a dot (.) and any blanks stripped off
      mr = irlong(croot,40)

      md = 0
      if(lpath) then
         read(pth_unit,'(a40)',err = 30) ctemp
         md = iclong(ctemp,40)
      end if
30    continue
      if(md.gt.0) then
         cfin = ctemp(1:md)//'/'//croot(1:mi)
      else
         cfin = croot(1:mi)
      end if

      if(l_asc) then
         open(unit=in_unit,file=cfin,form='formatted',status='old')
         if(.not.l_clk_hd) then
ccc         need to open a clock reset file ... set name here
            if(md.gt.0) then
               cfclk = ctemp(1:md)//'/'//croot(1:mr)//'.clk'
            else
               cfclk = croot(1:mr)//'.clk'
            endif
         endif
ccc      asc_rec will read clock reset info, find nch, position file at
ccc      begining of data
         call asc_rec(nch,time_int) 
      else
ccc      open binary file and read header (to find # of data channels)
         open(unit=in_unit,file=cfin,form='unformatted',access='direct',
     &   recl=256,status='old')
         call rdhd(nch)
         close(in_unit)

ccc      reopen binary file with correct record length
c        integer*4 and integer*2 binary input files are supprted. 
c        BYTES switches between
c        the two modes and the record length is computed using BYTES.
c        bytes is set equal to  parameter NBYTES (which is defined 
c        in PARAMS.INC) in dnff.f unless it isn't changed with the 
c        command line option -b ; See also rbblk
         irecl = (bytes*nblk+8)*nch+8
         open(unit=in_unit,file=cfin,form='unformatted',access='direct',
     &      recl=irecl,status='old')
      endif

c        system parameter directory (In all cases blank field => localdirectory)

      md = 0
      if(lpath) then
         read(pth_unit,'(a40)',err = 40) ctemp
         md = iclong(ctemp,40)
      end if
40    continue
      if(ctemp(1:8).eq.'standard') then
         cfsp='standard.sp'
      else
         if(md.gt.0) then
            cfsp = ctemp(1:md)//'/'//croot(1:mr)//'.sp'
         else
            cfsp = croot(1:mr)//'.sp'
         end if
      end if

      md = 0
      if(lpath) then
         read(pth_unit,'(a40)',err = 50) ctemp
         md = iclong(ctemp,40)
      end if
50    continue
      if(md.gt.0) then
         cfbr = ctemp(1:md)//'/'//croot(1:mr)//'.bad'
      else
         cfbr = croot(1:mr)//'.bad'
      end if

       
c       decset file (full path name)
      md = 0
      if(lpath) then
         read(pth_unit,'(a80)',err = 60) cfdecset
         md = iclong(cfdecset,80)
      end if
60    continue
      if(md.eq.0) then
         cfdecset = 'decset.cfg'
      end if

c       pwset file (full path name)
      md = 0
      if(lpath) then
         read(pth_unit,'(a80)',err = 70) cfpwset
         md = iclong(cfpwset,80)
      end if
70    continue
      if(md.eq.0) then
         cfpwset = 'pwset.cfg'
      end if

c       output directory
      md = 0
      if(lpath) then
         read(pth_unit,'(a40)',err = 80) ctemp
         md = iclong(ctemp,40)
      end if
80    continue
        if (nch.lt.10) then
           write (cfn,'(a1,i1)') 'f',nch
           lcfn = 2
        else if (nch.ge.10) then
           write(cfn,'(a1,i2)') 'f',nch
           lcfn = 3
        endif
      if(md.gt.0) then
         cfout = ctemp(1:md)//'/'//croot(1:mr)//'.'//cfn(1:lcfn)
      else
         cfout = croot(1:mr)//'.'//cfn(1:lcfn)
      end if
      close(pth_unit)
      return
      end
c______________________________________________________________________
c
      subroutine rdhd(nch)
      include 'iounits.inc'
      include 'input.inc'

      character*256 c
      character*80 comment
      common /FORMIN/irec

c    header read routine for binary 2-byte integer cleaned data
c     files
      irec = 1
      read(in_unit,rec = irec) c
      read(c,100) cfile,nch,dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment
100   format(5x,a80,5x,i2,4x,f12.4,4x,i2,4x,5i2,4x,i10,a80)

ccc   THIS IS PROBABLY WRONG:
ccc      time_int = iclk0*dr
ccc   I THINK IT SHOULD BE
      time_int = 0
ccc    (ACTUALLY: we need to think about/clarify the meaning of iclk0 in
ccc       binary data files!!!!!)
      sampfreq = 1./dr
      return
      end
c______________________________________________________________________
c
      subroutine rdblk(npts,ix,ierr,nch)

      include 'iounits.inc'
      include 'input.inc'
      include '../include/datsz.inc'
      parameter(nmx = nblk*nchmx)

cnew      character*80 cfile
cnew      logical lfirst,lfmiss
cnew      real sampfreq
cnew      common /inblk/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
cnew     &   ioff,inunit,nmsmx,lfirst,itf,lfmiss,cfile,time_int
cnew     &     sampfreq

c     BINARY TWO BYTE INTEGER VERSION

      integer*4 ix4(nmx)
      integer*2 ix2(nmx)
      integer ix(0:nch,*),it1,npts,ix0(nchmx),msval,ioff,ierr  !ioin,
     &   ,iscl(nchmx)
      common /FORMIN/irec

      n = nch*nblk
      irec = irec + 1

c     depending on BYTES read into an integer*2 or integer*4 array
      if (bytes .eq.2) then
         read(in_unit,rec=irec,err=99) it1,nt,(ix0(j),iscl(j),j=1,nch),
     &        (ix2(k),k=1,n)
      elseif (bytes .eq .4) then
         read(in_unit,rec=irec,err=99) it1,nt,(ix0(j),iscl(j),j=1,nch),
     &        (ix4(k),k=1,n)
      endif

      if(lfirst) then
         nmiss = 0
      else
         nmiss = it1 - itf
         if(nmiss.lt.0) then
            ierr = -2
            return
         end if
         if(nmiss .gt. nmsmx) then
            ierr = -3
            nmiss = 0
         end if
      end if

      do i = 1,nmiss
      ix(0,i) = itf - 1 + i + ioff
         do j = 1,nch
            ix(j,i) = msval
         enddo
      enddo

      npts = nt + nmiss
      do i = 1,nt
         ix(0,i+nmiss) = i -1 + it1 + ioff
         jj = nch*(i-1)
         do j = 1,nch
            if (bytes .eq .2) then
               ix(j,i+nmiss) = ix0(j) + iscl(j)*ix2(jj+j)
            elseif (bytes .eq. 4) then
               ix(j,i+nmiss) = ix0(j) + iscl(j)*ix4(jj+j)
            endif
         enddo
      enddo

      itf = ix(0,npts) + 1 - ioff
      ierr = 0

      return                              

99    continue
      ierr = -1
      return
      end
c___________________________________________________________
c
      function iclong(croot,n)
ccc   find positon of last non-blank characters in a string
ccc    (i.e., end of string, hopefully)
      character*1 croot(n)
      integer iclong
      do 10 i = n,1,-1
      if(ichar(croot(i)).ne.32) then
         iclong = i
         return
      end if
10    continue
      iclong = 0
      return
      end

c___________________________________________________________
c
      function icfirst(croot,n)
ccc   find number of blank characters at begining of a string
      integer icfirst
      character*1 croot(n)
      do 10 i = 1,n
      if(ichar(croot(i)).ne.32) then
         icfirst = i
         return
      end if
10    continue
      icfirst = n
      return
      end
c___________________________________________________________
c
      function irlong(cstr,n)
ccc   find number of characters before first '.' or ' ' in a string
      integer irlong
      character*1 cstr(n)
      do 10 i = 1,n
      if((cstr(i).eq.' ').or.(cstr(i).eq.'.')) then
         irlong = i-1
         return
      end if
10    continue
      irlong = n
      return
      end
