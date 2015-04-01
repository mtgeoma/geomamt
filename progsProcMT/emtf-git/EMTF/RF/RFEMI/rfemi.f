      program rfemi
c     program to reformat EMI files into standard binary direct access

c*****************************************************************************
ccc########################################################################
ccc   change for number of channels
c     number of channels = nch
      include '../../include/nchmx.inc'
      parameter (nblkmx=5000)
      character*10 chid(nchmx)
      
ccc########################################################################

c******parameters to change for output file format
      parameter(nchmxp1 = nchmx+1, ntotmx = 1000)
        
      integer ix(nchmxp1*ntotmx)
      character*80 comment,cfile,ctemp
      integer iyr,imo,iday,ihr,imin,isec,iclk0
      real dr
      common /HEADER/dr,iyr,imo,iday,ihr,imin,isec,iclk0,comment
      parameter(nsegmx = 100)
      integer nseg,iseg(2,nsegmx),ixsv(nchmx,nblkmx)
      real xb(nchmx),xsc(nchmx)
      logical lend,lfirst,lsp,lhd
      character*1 ans
      character*3 csta
      character*80 cfout

      lhd = .false.
      narg = iargc()
      do k = 1,narg
         call getarg(k,ctemp)
         if(ctemp(1:2).eq.'-h') lhd = .true.
      enddo

1     continue
      print*,'station id'
      read(5,'(a3)') csta
      write(0,*) csta
      ii = 0
      nseg = 0

      print*,'enter output file name'
      read(5,'(a80)') cfout
 
c>>>>>>>>>>>   initialize input
      inunit = 3
      lsp = .true.
      call ininit(inunit,cfout,lsp,lhd,nch,chid)
      do i = 1,nch
         chid(i) = chid(i)(1:2)//':'//csta
      enddo
      nchp1 = nch+1

c>>>>>>>>>>>   initialize output
      iout = 1
      call outinit(iout,nch,nt,cfout)

c      print*,'nt = ',nt
3      continue 

      call frstdat(inunit)
      nseg = nseg + 1
      lfirst = .true.
c>>>>>>>>>>  main loop starts here; get ntot records
40      continue  
c       get a block of data
      j0 = 0
      lend = .false.
      call indo(inunit,ix,nt,ngot,j0,nchp1,lend)
      if(lfirst) then
         it1 = ix(1)
         iseg(1,nseg) = it1
         ir = it1
      end if
      lfirst = .false.
c      save one output record from each block to estimate scales, mean
           ii = ii + 1
           do 150 k = 1,nch
150           ixsv(k,ii) = ix(k+1)

c            see if a new segment should be initialized
           do i0 = 1,ngot-nt+1,nt 
              call outchk(ix(1+nchp1*(i0-1)),nch,nt,iout,ir,iseg,nseg)
              i1 = i0
           enddo
           i0 = i1 + nt - 1
           j0 = ngot - i0
           do i = i0+1,ngot
              do k = 1,nchp1
                 ix(k+nchp1*(i-i0-1)) = ix(k+nchp1*(i-1))
              enddo
           enddo
        if(.not.lend) then
           go to 40
        else
           iseg(2,nseg) = ir-1
        end if 

c     append another data segment to current output file?
      print*,'another input file? (append to current output file)'
      read(5,'(a1)') ans
      write(0,*) 'ans = ',ans
      if( (ans.eq.'y') .or. (ans.eq.'Y') ) then
         lsp = .false.
         call ininit(inunit,cfile,lsp,nch,chid)
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
      write(0,*) 'ans = ',ans
      if(ans.eq.'y'.or.ans.eq.'Y') go to 1
      stop
      end
