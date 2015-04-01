c______________________________________________________________________
c
      subroutine cininit(inu,nch)
      include '../../include/datsz.inc'
      character*80 cfile
      character*80 comment
c       cfsp, cfbr,cfdecset,cfpwset,cdirout are path/file names for
c       system params, bad recs,decset,pwset, and output file
      common /inblk/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
     &  ioff,inunit,itf,comment

      ioff = 0

      inunit = inu
      print*,'enter input file name'
      read(5,'(a40)') cfile

      open(unit=inunit,file=cfile,form='unformatted',access='direct',
     &   recl=256)
      call rdhd(nch)
      write(0,*) 'nch = ',nch
      close(inunit)
      irecl = (2*nt0+8)*nch+8
      open(unit=inunit,file=cfile,form='unformatted',access='direct',
     &   recl=irecl)
      return
      end
c______________________________________________________________________
c
      subroutine rdhd(nch)

      integer iyr,imo,iday,ihr,imin,isec,iclk0
      character*256 c
      character*80 cfile                           
      character*80 comment
      common /FORMIN/irec
      common /inblk/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
     &  ioff,inunit,itf,comment

c    header read routine for binary 2-byte integer cleaned data
c     files
      irec = 1
      read(inunit,rec = irec) c
      read(c,100) cfile,nch,dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment
100   format(5x,a80,5x,i2,4x,f12.4,4x,i2,4x,5i2,4x,i10,a80)
      return
      end
c______________________________________________________________________
c
      subroutine rdblk(nt,nch,it1,ix,ierr)

      include '../../include/datsz.inc'
      parameter(nmx=nchmx*nt0)

      character*80 comment
      common /inblk/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
     &   ioff,inunit,itf,comment

c     BINARY TWO BYTE INTEGER VERSION

      integer*2 ix2(nmx)
      integer ix(nch,*),it1,itf,ix0(nchmx),ioff,ierr
     &   ,iscl(nchmx)
      common /FORMIN/irec

      n = nch*nt0
      irec = irec + 1
      read(inunit,rec=irec,iostat=ier) it1,nt,(ix0(j),iscl(j),j=1,nch),
     &   (ix2(k),k=1,n)
      if(ier .lt. 0) go to 99

      do i = 1,nt
         jj = nch*(i-1)
         do j = 1,nch
            ix(j,i) = ix0(j) + iscl(j)*ix2(jj+j)
         enddo
      enddo
      ierr = 0
      return                              

99    continue
c       print*,'error reading input file'
      ierr = -1
      return
      end
