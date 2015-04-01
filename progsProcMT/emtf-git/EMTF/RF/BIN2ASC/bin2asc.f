      parameter (nchmx = 11)
      character*80 comment
      character*1 ans
      integer nch,iyr,imo,iday,ihr,imin,isec,iclk0,ioff,
     & inunit,itf
      integer ix(nchmx,1000),it1
      common /inblk/dr,iyr,imo,iday,ihr,imin,isec,iclk0,
     &  ioff,inunit,itf,comment
      logical lsamp 
      character*80 arg

      lsamp = .false.
      maxblks = 1000000
      narg = iargc()
      do k = 1,narg
         call getarg(k,arg)
         if(arg(1:2).eq.'-M') then
ccc         change maximum number of blocks to read
ccc         use this to format a portion of the file
            read(arg(3:80),*) maxblks    
         elseif(arg(1:2).eq.'-S') then
ccc         output sample numbers in ASCII data file
            lsamp = .true.
         endif
      enddo
      inu = 2
      iout = 3
10    continue
      call cininit(inu,nch)
      call outinit(iout,nch)

      do m = 1,maxblks
         call rdblk(nt,nch,it1,ix,ierr)
         if(ierr.lt.0) then 
            close(inu)
            close(iout)
            print*,'convert another file?'
            read(5,'(a1)') ans
            if((ans.eq.'y').or.(ans.eq.'Y')) go to 10
            stop
         else
            call outdo(it1,nt,nch,ix,lsamp,iout)
         end if
      enddo
      end
