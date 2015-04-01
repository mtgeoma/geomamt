c______________________________________________________________________
c
      subroutine outinit(iout,nch)

      character*80 cfout

c  commonblock FORMOUT contains info about output format
c    (common to outinit and outdo)
      character*40 cf,cfhd
      common /FORMOUT/cf,cfhd

c  common block HEADER contains info passed to output file, but not
c    used in main program
      character*80 comment
      integer iyr,imo,iday,ihr,imin,isec,iclk0,ioff,inunit,itf
      real dr
      common /inblk/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
     &  ioff,inunit,itf,comment

      print*,'enter output file name'
      read(5,'(a80)') cfout
      open(unit=iout,file=cfout)

      return
      end
c______________________________________________________________________
c
      subroutine outdo(it1,nt,nch,ix,lsamp,iout)

c  commonblock FORMOUT contains info about output format
c    (common to outinit and outdo)
      character*40 cf,cfhd
      character*80 comment
      common /FORMOUT/cf,cfhd
      common /inblk/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
     &  ioff,inunit,itf,comment
      logical lsamp

      parameter (ntmx=1000)
      integer ix(nch,*),it1,nt

      if(lsamp) then
         do j =  1,nt
            write(iout,'(i10,10i7)') it1+j-1,(ix(i,j),i=1,nch)
         enddo
      else
         do j =  1,nt
            write(iout,'(10i7)') (ix(i,j),i=1,nch)
         enddo
      endif
 
      return
      end
c______________________________________________________________________
c
      subroutine wrhdbin(iout,cfile,nch,dr,iyr,imo,iday,ihr,imin,
     & isec,iclk0,irecl,comment,nseg,iseg,xb,xsc,chid)

ccc   here this is just a dummy routine ... nothing happens when this is
ccc   called ... could easily be changed to have some sort of header output
ccc   if desired

      return
      end
