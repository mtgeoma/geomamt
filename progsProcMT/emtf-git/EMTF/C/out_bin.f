      subroutine outinit(iout,nch,nt)

c   THIS VERISON IS FOR BINARY (2-byte integer)
c         DIRECT ACCESS OUTPUT FILES:
C    ONE DIRECT ACCESS RECORD IS A BLOCK IN THE PAKED ASCII
C    CHARACTER FORMAT

      include '../include/datsz.inc'
      parameter(nmx=nchmx*nt0)

      character*80 cfout

c  commonblock FORMOUT contains info about output format
c    (common to outinit and outdo)
      common /FORMOUT/irec,ix2

c  common block HEADER contains info passed to output file, but not
c    used in main program
      character*80 comment
      integer iyr,imo,iday,ihr,imin,isec,iclk0
      integer*2 ix2(nmx)
      real dr
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment

      irecl = (2*nblk+8)*nch+8
      nt = nblk

      print*,'enter output file name'
      read(5,'(a80)') cfout
      open(unit=iout,file=cfout,form='unformatted',access='direct',
     &   recl=irecl)

      print*,'enter header (80 character max)'
      read(5,'(a80)') comment

      irec = 1

      return
      end
c______________________________________________________________________
c
      subroutine wrhdbin(iout,cfile,nch,dr,iyr,imo,iday,ihr,imin,
     & isec,iclk0,irecl,comment,nseg,iseg,xb,xsc,chid)

c     like wrhd, but for binary file format

      character*256 c
      character*80 cfile
      character*80 comment

      character*10 chid(nch)
      integer nseg,iseg(2,nseg)
      real xb(nch),xsc(nch)

      write(0,*) 'in wrhdbin, nch = ',nch
 
      do 5 i = 1,256
5     c(i:i) = ' '
      write(c,102) cfile,nch,dr,iyr,imo,iday,ihr,imin,isec,iclk0
     &,comment

      write(iout,rec=1) c,nseg,(iseg(1,k),iseg(2,k),k=1,nseg),
     &  (chid(l),xb(l),xsc(l),l=1,nch)

102   format('fn=  ',a80,',nch=',i2,',dr=',f12.4,',yr=',
     & i2,',t0=',5i2,',i0=',i10,a80)

      return
      end
c_____________________________________________________________________
c
      subroutine outdo(it1,nt,nchp1,ix,iout)

      include '../include/datsz.inc'
      parameter(nmx=nchmx*nt0)


c  commonblock FORMOUT contains info about output format
c    (common to outinit and outdo)
      common /FORMOUT/irec,ix2

c       as for wrblkvg, but output is in binary integer*2 fixed record length
c       blocks

      integer ix(nchp1,*),it1,nt,ixmin(100),ixmax(100),iscl(100)
      integer*2 ix2(nmx)
      logical lix1

c         find minimum, maximum for each channel
      do 2 i = 2,nchp1
      ixmin(i) = ix(i,1)
      ixmax(i) = ix(i,1)
         do 1 j = 2,nt
         ixmax(i) = max(ixmax(i),ix(i,j))
1        ixmin(i) = min(ixmin(i),ix(i,j))
      irng = ixmax(i)-ixmin(i)
      if(irng.le.32768) then
         iscl(i) = 1
      else
         iscl(i) = nint(irng/32768.+.5)
      end if
2     continue
 
c          subtract min, rescale, output
      ii = 0
      do 10 i = 1,nt
      do 10 j = 2,nchp1
      ii = ii + 1
      ix2(ii) = nint(float(ix(j,i) - ixmin(j))/iscl(j))
10    continue
      irec = irec + 1
      write(iout,rec=irec) it1,nt,(ixmin(j),iscl(j),j=2,nchp1),
     &      (ix2(i),i=1,ii)

      return
      end
