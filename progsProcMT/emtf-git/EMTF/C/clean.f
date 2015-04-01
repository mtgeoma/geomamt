         program clean
c       program to do a preliminary cleanup of data points (e.g. remove
c       parity errors)

c**********************************************************************
c     number of channels = nch
c       NEED TO CHANGE THIS TO RUN WITH DIFFERENT NUMBERS OF CHANNELS
      parameter (nch = 5)
c       NEED TO CHANGE THIS TO ALLOW FOR VERY LONG DATA FILES
c        nblkmx is maximum number of blocks
c           nblk (defined in ../datsz.h) is number of samples per block
c        (total number of samples in output file should be <= nblkmx*nblk)
      parameter (nblkmx=5000)
c**********************************************************************

c    n is number of points used to determine reasonable scale for data
      parameter (n=7)
      parameter(nchp1 = nch+1, ntotmx = 1000,
     1   nmx = n+n/2-1,nmx2= 2*nmx,nxplt=nch*2+1,nxplt1=nxplt-1)
        
      integer ix(nchp1,ntotmx),ixcln(nchp1,nmx2),itemp(nmx2)
     1,irecs(2),nmad(6),ixt(nchp1),ixplt(nxplt,ntotmx),ixsv(nch,nblkmx)
        
      parameter(nsegmx = 100)
      character*10 chid(5)
      integer nseg,iseg(2,nsegmx)
      real xb(nch),xsc(nch)

      real scale(6),pscale(6)
      logical lend,lplot,lclean
      character*80 cftemp
      character*7 stname
      character*80 cfplt
      character*1 ans,c1
      character*4 c4
      character*24 cphd
      character*2 cx(nchp1,ntotmx)
      character*3 csta
      
      data scale/1.,2.,2.,2.,2.,2./
      data pscale/1.,1.,1.,1.,1.,1./
c      data nmad/3,8,8,8,10,10/
      data nmad/3,800,800,800,1000,1000/
      data chid/'H','D','Z','E','N'/
ccc   lclean is true to clean data with median filtering
ccc   if false, just reformat and output
      data lclean/.true./

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c   BEGIN EXECUTABLE

1     continue
      print*,'station id'
      read(5,'(a3)') csta
      do 2 i = 1,nchmx
2        chid(i) = chid(i)(1:1)//':'//csta
      ir = -1
      nseg = 0
      ii = 0

c>>>>>>>>>>>   initialize input
      inunit = 2
      call ininit(inunit)
c>>>>>>>>>>>   initialize output
      iout = 1
      call outinit(iout,nch,nt)
      ntot = nt+n-1
c>>>>>>>>>>>   initialize plotting
c      print*,'enter plot file name'
c      read(5,'(a80)') cfplt
c      open(unit=3,file=cfplt)

3      continue 
c>>>>>>>>read in initial records and pad left end of time series
      lend = .false.
      n2 = n/2
      nmed = n2+1
      call indo(inunit,ix,n,ngot,1,nchp1,lend)
      do 5 i = 1,n
         itemp(i) = ix(1,i)
5        continue                                                         
c   pad record numbers
      call isort(itemp,n)
c   irec is first (cleaned) record number
      irec = itemp(nmed) - n2
      do 7 i = 1,n2
7        ix(1,i) = irec - n2 + i -1

c      now pad other channels
       do 18 j = 2,nchp1
          do 15 i = 1,n
15           itemp(i) = ix(j,i)
          call isort(itemp,n)
          do 16 i = 1,n2
16           ix(j,i) = itemp(nmed)
18        continue

c>>>>>>>>>>> done with initial record input and padding
c      reposition input file at start of first record
c     and re-initialize pointers

      call frstdat(inunit)
      ipoint = nmed
      i1 = nmed
      lrec = irec - n
      nrec = irec + nt

      iyiyi = 0
c>>>>>>>>>>  main loop starts here; get ntot records for cleaning
40      continue  
       iyiyi = iyiyi + 1
      call indo(inunit,ix,ntot,ngot,ipoint,nchp1,lend)
      i2 = ngot - n -n2 + 1
       
45    iplt = 0
      if(lclean) then
      do 100 i = i1,i2,n
         iplt = iplt + 1
c       preliminary screening
         i1n = i + n - 1
         do 550 j=1,nchp1
            ibar = 0
            do 500 k = i,i1n
500            ibar = ibar + ix(j,k)
            ibar = ibar/n
            imax = abs(ix(j,i)-ibar)
            do 510 k = i+1, i1n
510            imax = max(imax,abs(ix(j,k)-ibar))
            ixp1 = (j-1)*2
            if(j.eq.1) ixp1 = 1
            if(imax.gt. nmad(j)) then
c          maximum deviation from mean is too large; check points more carefully\
c          and fix if necessary; note that cleaned data comes out in array ixcln
               call cleanrec(ix(1,i),n,nchp1,j,scale,ixcln,
     1                      ixplt(ixp1,iplt),pscale ,nmad(j),nmed)
c               move cleaned data back into array ix
                  do 520 k = 1,n
520                 ix(j,i-1+k) = ixcln(j,k) 
              else
                 ixplt(ixp1,iplt) = ibar
                 ixplt(ixp1+1,iplt) = imax*pscale(j)
              end if
550        continue
100     continue
        else
           if(i2.lt.i1) then
              iplt = 0
           else
              iplt = (i2-i1)/n + 1
           endif
        endif    ! lclean

        if(.not. lend) then
           it1 = ix(1,i1)
           npts = iplt*n
c      save one output record from each block to estimate scales, mean
           ii = ii + 1
           do 150 k = 1,nch
150           ixsv(k,ii) = ix(k+1,i1)

c            see if a new segment should be initialized
           if(ir.ne.it1) then
              print*,ir,nseg,it1,npts
              nseg = nseg + 1
              iseg(1,nseg) = it1
              if(nseg.ne.1) iseg(2,nseg-1) = ir-1
              ir = it1+npts
           else
              ir = ir + npts
           end if
           call outdo(it1,npts,nchp1,ix(1,i1),iout)
ccc        move left over data
           in2 = i + n2 + 1 - n
           do 120 l = in2,ngot
              ii = l - in2 + 1
              do 120 j = 1,nchp1
120           ix(j,ii) = ix(j,l)

c       reset pointers
           irecs(1) = irec
           irecs(2) = irec+nt-n
           ipoint = ii + 1
           irec = nrec
           lrec = nrec - n
           nrec = nrec + nt
           go to 40
        end if

200     continue
        if(lclean) then
c       finish up; output last few points
c       pad end of series before calling cleaning routine
c       first do time channel

        npoints = ngot - ( i - 1 )
        if(ngot .lt. n) go to 300
        do 205 k = 1,n
           i1 = ngot - k + 1
205        itemp(k) = ix(1,i1)
        call isort(itemp,n)
        irec0 = itemp(nmed) + n2
        do 210 k = 1,n2
210        ix(1,ngot + k) = irec0 + k

           do 250 j = 2,nchp1
              do 215 k = 1,n
              i1 = ngot -k + 1
215           itemp(k) = ix(j,i1)
           call isort(itemp,n)
           do 220 k = 1,n2
220           ix(j,ngot+k) = itemp(nmed)
250        continue

        i1 = ngot - n + 1
           do 265 j = 1,nchp1
              ixp1 = (j-1)*2
              if(j.eq.0) ixp1 = 1
              call cleanrec(ix(1,i1),n,nchp1,j,scale,ixcln,
     .             ixplt(ixp1,1),pscale,nmad(j),nmed)
               do 260 k = 1,n
260               ix(j,i1-1+k) = ixcln(j,k)
265           continue

        if(npoints .gt. n) then
           i1 = i    
           do 275 j = 1,nchp1
              ixp1 = (j-1)*2
              if(j.eq.0) ixp1 = 1
              call cleanrec(ix(1,i1),n,nchp1,j,scale,ixcln,i
     .              xplt(ixp1,1),pscale,nmad(j),nmed)
              do 270 k = 1,n
270              ix(j,i1-1+k) = ixcln(j,k)
275           continue
        end if

c     write out last few data points
      endif    ! lclean
300     continue
      npts = ngot - n2
      it1 = ix(1,nmed)
      write(0,*) 'npts,ii',npts,ii
      do 310 k = 1,nch
310     ixsv(k,ii) = ix(k+1,i1)

c           see if a new segment should be initialized
        if(ir.ne.it1) then
           print*,ir,nseg,it1,npts
           nseg = nseg + 1
           iseg(1,nseg) = it1
           if(nseg.ne.1) iseg(2,nseg-1) = ir-1
           ir = it1+npts
        else
           ir = ir + npts
        end if

      call outdo(it1,npts,nchp1,ix(1,nmed),iout)
 
c        append another data segment to current output file?
      print*,'another input file? (append to current output file)'
      read(5,'(a1)') ans
      if( (ans.eq.'y') .or. (ans.eq.'y') ) then
         call ininit(inunit)
         go to 3
      else
         iseg(2,nseg) = ir-1
         do 420 k = 1,5
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
