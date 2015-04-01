      subroutine ininit(inunit,cfout,lsp,lhd,nch,chid)

c     initializes input file: reads header, passes header data
c     to ouput initialization, computes clock zero etc;

      include '../../include/nchmx.inc'
      include '../../include/four_byte.inc'
ccc   NOTES:

ccc   if lsp = .true. a system parameter file is made

ccc   This routine expects to find a "clock reset" file, which tells
ccc   the time of the first 

      character*80 cfile,cfhd,cfout,ctemp
      character*40 cformat
      character*20 csp
      character*2 chid

c  common block HEADER contains info passed to output file, but not
c    used in main program
      character*80 comment
      integer iyr,imo,iday,ihr,imin,isec,iclk0,iord(nchmx)
      real dr,clock(2)
      logical lsp
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
     &   comment

      integer mday(12)
      data mday/31,28,31,30,31,30,31,31,30,31,30,31/
      integer ir0,c_in_size
      common /FORMBLK/ir0,irelrc,irec0,nt,nbl,iord
      common /CBLK/irec,cfile,nskip

      parameter (c_in_size = 2000)
      character*1 chead(81,100),ct,c_in(c_in_size)
      logical lcr,lhd
      integer ich(100)
      data icr/13/
      data ilf/10/

1     continue
      print*,'input file name'
      if(l_4byte) then
         irecl = c_in_size/4
      else
         irecl = c_in_size
      endif
      read(5,'(a80)') cfile
      open(unit=inunit,file=cfile,form='unformatted',access='direct',
     &    recl=irecl,status='old')

ccc   read in a bunch of characters, including full EMI header
      read(inunit,rec=1) c_in
ccc   then sort out the contents ...
      ir = 0
ccc   read header as a series of ASCII strings; sort out contents
ccc   (if desired) later in routine mksp     
      ihd = 1
      ich(ihd) = 0
      lcr = .false.
10    continue
      if(ich(ihd).eq.81) go to 20
         ir = ir+1
ccc      now copy one character at a time into ct, just as if we were reading the
ccc      file one character at a time   .... that's the qd fix
         ct = c_in(ir)
c         read(inunit,rec=ir) ct
         if(lcr) then
            if(ichar(ct).eq.ilf) then
ccc             end of header record ihd
c               print*,'ihd,ich(ihd) = ',ihd,ich(ihd)
               ihd = ihd + 1
               ich(ihd) = 0
               lcr = .false.
            else
               ich(ihd) = ich(ihd) + 2
               chead(ich(ihd),ihd) = ct
               lcr = .false.
               write(0,*) 'CR not followed by LF'
            end if
            go to 10
         else if(ichar(ct).eq.icr) then
            lcr = .true.
            go to 10
         else
            lcr = .false.
            ich(ihd) = ich(ihd) + 1
            chead(ich(ihd),ihd) = ct
            go to 10
         end if
20    continue
c        ir0 marks last character in header block;
c         first data point is in bytes ir0+1, ir0+2
      print*,'ir = ',ir
      ir0 = ir - 81
      nskip = ir0
      nhd = ihd-1
      if(nhd.eq.0) then
         print*,'ERROR: stopping'
         stop
      endif

c     print out header info
      if(lhd) then
ccc      output EMI file header verbatim
         ctemp = cfout
         ctemp(1:2) = 'hd' 
         write(0,*) 'ctemp = ',ctemp
         open(unit=77,file=ctemp,status='unknown')
         do i = 1,nhd
            write(77,'(80a1)') (chead(j,i),j=1,ich(i))
         enddo
         close(77)
      endif
 
c>>>>>>>>>>>>>>>>    read header info from clock reset file
      print*,'enter clock reset file name'
      read(5,'(a80)') cfhd
      open(unit = 10, file = cfhd,status='unknown')
ccc   skip sampling rate line (so same file can be used for different
ccc   sampling bands
      read(10,*)

c*instrument clock zero time - year month, day, time (ut) [ hour,min,sec]
      read(10,*) iyr,imo,iday,ihr,imin
cc*universal clock zero time - year month, day, time (ut) [ hour,min,sec]
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

      if(lsp) call mksp(cfout,chead,nhd,ich,isfreq,nch,clock,chid,dr)

c     compute julian day number for instrument clock zero
        jday = 0
        do 120 i = 1,imo-1
        if( (mod(iyr,4).eq.0) .and. (i.eq.2) ) then
           jday = jday+mday(i) + 1
        else
           jday = jday + mday(i)
        end if
120     continue
        jday = jday + iday

cc     compute julian day number for universal clock zero
        jdayu = 0
        do 125 i = 1,imou-1
        if( (mod(iyru,4).eq.0) .and. (i.eq.2) ) then
           jdayu = jdayu+mday(i) + 1
        else
           jdayu = jdayu + mday(i)
        end if
125     continue
        jdayu = jdayu + idayu

c   record number (relative to universal clock zero) of instrument
c      clock zero

      jday = jday-jdayu + 1
      ihr = ihr - ihru
      imin = imin - iminu
      isec = isec - isecu
c         iclk0 is absolute record number of external clock zero
      iclk0 = (((((jday-1)*24+ihr)*60+imin)*60+isec))*isfreq

ccc   now get starting minute (external clock) of time series from first
ccc   block in data file

      ir = ir0 + nch*(nbl*2+2)
      irec1 = int((ir+1)/c_in_size)+1
      irec2 = int((ir+11)/c_in_size)+1
      if(irec1.ne.irec2) then
           write(0,*) 'Timing will be wrong ...'
           write(0,*) 'Look at source code input.new.emi...'
ccc         If you get this message it is because
ccc         we are reading in characters from the beginig
ccc         of the EMI file in chunks, and then searching 
ccc         through character strings looking for a particular
ccc         pattern ... in this case the 'OK' at the end of
ccc         the first data block so we can read off the
ccc         external clock time... irec1 .ne. irec2 means
ccc         that with the chosen block size for reading the
ccc         charactrers we want to search span two blocks,
ccc          and making the program general enough to deal with
ccc         this seems ridiculous, considering that everything
ccc         (including the EMI file format) is such a kluge.
      endif
      read(inunit,rec=irec1) c_in
      close(inunit)
      ir = ir - c_in_size*(irec1-1)
      do 205 i = 1,11
         csp(i:i) = c_in(ir+i)
         if(csp(i:i).eq.'O') go to 210
      print*,'i,ir,csp(i:i)',i,ir,csp(i:i)
205      continue
c       il is length of integer string giving external clock start time
c             for data file
c          il + 2 is length of clock record at end of block
210   il = i-1
      if(il.lt.10) write(cformat,211) il
      print*,'il = ',il
   
211   format('(i',i1,')')

c     nt is total number of bytes in each data block
      nt  = nch*(nbl*2+2)+il+4
      if((il .gt. 0).and.(il.lt.10)) then
         read(csp,cformat) istart
         print*,'istart = ',istart
      else
         istart = 0
      endif

c       irec0 is universal record number of first sample in data file
c       first, assuming clocks are OK
      print*,'iclk0 = ',iclk0
      irec0 = iclk0 + istart*60*isfreq
      print*,'istart,isfreq,istart*60*isfreq',istart,isfreq,
     &    istart*60*isfreq
      irec0 = irec0 + iroff
      print*,'irec0 = ',irec0

      return
      end
c_____________________________________________________________________c
      subroutine frstdat(inunit)
   
CCC   version to initialize for C reading routine
ccc   close inunit, 
c        position input file at firstheader record

      include '../../include/nchmx.inc'
      character*80 comment,cfile,cformat
      integer iyr,imo,iday,ihr,imin,isec,iclk0,iord(nchmx)
      integer*2 file_open_c
      real dr
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment
      common /CBLK/irec,cfile,nskip
      common /FORMBLK/ir0,irelrc,irec0,nt,nbl,iord

      close(inunit)
      n_cfile = 80
      do while(cfile(n_cfile:n_cfile) .eq. ' ')
         n_cfile = n_cfile - 1
      enddo
      open(unit=11,file = 'file_name_temp',status='unknown')
      if(n_cfile.lt.10) then
         write(cformat,11) n_cfile
11       format('(a',i1,')')
      else
         write(cformat,12) n_cfile
12       format('(a',i2,')')
      endif
      write(11,cformat) cfile(1:n_cfile)
      close(11)
      ierr = file_open_c(nskip)

      if(ierr .ne. 0) then
         write(0,*) 'Error opening file for C input'
         write(0,*) 'ierr = ',ierr
      endif
      irec = irec0
      return
      end
c_____________________________________________________________________c
      subroutine indo(inunit,ix1,n,ngot,j0,nchp1,lend)

ccc    this version calles C routine rdblk to get a block of EMI
ccc     MT-1 data and position file at start of next block
c      inunit = input unit number
c      ix1 = data array (integer)
c      n = number of points to try returning in ix1
c      ngot = number of points actually returned
c      nch = number of data channels
c      j0+1 = starting index in array ix1 (smaller
c         indices in array already have data)
c      lend = .true. if EOF is reached (logical variable)
c      
      include '../../include/datsz.inc'
      parameter(nmx=nchmx*nt0)
      parameter(nccmx = nmx + 20)

      integer ix1(nchp1,*),rdblk,j0
      integer*2 ix(nmx)
      character*80 cfile
      character*1 cc(2,nccmx)
      logical lend

ccc   iord(nch) gives order to store data channels in: put first channel in file
      integer iord(nchmx)
ccc   irec gives next record to read
      integer irec
      common /CBLK/irec,cfile,nskip
      common /FORMBLK/ir0,irelrc,irec0,nt,nbl,iord
        
      lend = .false.
      nch = nchp1-1
      ierr = rdblk(nch,nbl,ix)
      if(ierr .eq. 0) then
         ngot = 0
         lend = .true.
         write(0,*) 'End of File'
         return
      else if(ierr .eq. -1) then
         write(0,*) 'Gap'
         irec = irec + 2*nbl
      endif
      do j = 1,nbl
         irec = irec + 1
         ix1(1,j) = irec
         do i = 1,nch
            ix1(iord(i),j0+j) = ix(i+nch*(j-1))
         enddo
         ngot = nbl+j0
      enddo
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
