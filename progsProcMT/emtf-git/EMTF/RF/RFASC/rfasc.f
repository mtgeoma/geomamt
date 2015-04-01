      program rfemi
c     program to reformat data files into (my) standard binary direct access

c*****************************************************************************
      include 'nch.inc'
      integer ix(nchp1,ntotmx)

      character*80 comment,cfile,cfout
      integer iyr,imo,iday,ihr,imin,isec,iclk0
      real dr
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment

      parameter(nsegmx = 100)
      integer nseg,iseg(2,nsegmx),ixsv(nch,nblkmx)

      real xb(nch),xsc(nch)

      logical lend,lsp
      character*1 ans
      character*3 csta
      
1     continue
      print*,'station id'
      read(5,'(a3)') csta
      do 2 i = 1,nch
2        chid(i) = chid(i)(1:1)//':'//csta
c>>>>>>>>>>>   initialize output
      iout = 1
      print*,'enter output file name'
      read(5,'(a80)') cfout
      call outinit(iout,nch,nt,cfout)
c      write(0,*) 'nt,nch = ',nt,nch

c>>>>>>>>>>>   initialize input
      inunit = 3
      lsp = .false.
      call ininit(inunit,cfile,lsp)
      print*,'nt = ',nt
3      continue 

      call frstdat(inunit)
      inext = -1
      ii = 0
      nseg = 0
      ipoint = 1
c>>>>>>>>>>  main loop starts here; get ntot records
40      continue  
c       get a block of data
      lend = .false.
      call indo(inunit,ix,nt,ngot,ipoint,nchp1,lend)
c      write(0,*) 'ngot = ',ngot
      if(.not.lend) then
         it1 = ix(1,1)
c         write(0,*) 'it1 = ',it1
         if(it1.ne.inext) then
ccc         begin new segment
            if(nseg.gt.0) then
               iseg(2,nseg) = inext-1
            endif
            nseg = nseg + 1
            iseg(1,nseg) = it1
            inext = it1 + ngot
c            write(0,*) 'nseg,iseg,inext',nseg,iseg(1,nseg),inext
         else
            inext = inext + ngot
         end if
c      save one output record from each block to estimate scales, mean
         ii = ii + 1
         do k = 1,nch
            ixsv(k,ii) = ix(k+1,1)
         enddo
         call outdo(it1,ngot,nchp1,ix,iout)
         go to 40
      else
         if(ngot.gt.0) then
c            write(0,*) 'ngot = ',ngot
            it1 = ix(1,1)
            inext = inext + ngot
             
            call outdo(it1,ngot,nchp1,ix,iout)
         end if
         iseg(2,nseg) = inext -1
      end if 
c      write(0,*) 'ii = ',ii


c        append another data segment to current output file?
      print*,'another input file? (append to current output file)'
      read(5,'(a1)') ans
      if( (ans.eq.'y') .or. (ans.eq.'y') ) then
         lsp = .false.
         call ininit(inunit,cfile,lsp)
         go to 3
      else
         do 420 k = 1,nch
            xb(k) = 0
            do 410 i = 1,ii
410            xb(k) = xb(k)+ixsv(k,i)
            xb(k) = xb(k)/ii
            xsc(k) = 0.0
            do 415 i = 1,ii
415            xsc(k) = xsc(k) + abs(ixsv(k,i)-xb(k))
            xsc(k) = xsc(k)/ii
420         continue

         call wrhdbin(iout,cfile,nch,dr,iyr,imo,iday,ihr,imin,
     &     isec,iclk0,irecl,comment,nseg,iseg,xb,xsc,chid)
         close(iout)
      end if
      print*,'continue?'
      read(5,'(a1)') ans
      if(ans.eq.'y'.or.ans.eq.'Y') go to 1
      stop
      end
