      program dnff

      character*10 cdate
      character*5 cvers
      parameter (cvers = '5.1.1', cdate = '03/10/1998')

      include 'iounits.inc'
      include 'params1.inc'
      include 'decimate.inc'
      include 'input.inc'

      parameter (nt=nwmx/2,ntot= nwmx*nchmx,nchmx1=nchmx+1,
     &  nmax = 3*nwmx,nwsvmx = nwmx*4+15,nfmax=nwmx/2,nwmx1=nwmx+10
     &       , nchmx2 = 2*nchmx,nchmx21 = nchmx*2+1)
          
      real x(nchmx,nwmx1),samprate(ndmx),w(nwmx,ndmx)
     &    ,sc(nchmx),wsave(nwsvmx,ndmx),stcor(2),orient(2,nchmx),
     &    xx(ndmx,nchmx1,0:nmax),pspecl1(nchmx,nwmx,ndmx),
     &    areg(nchmx,10),rr(nwmx),ar(nwmx),bwo(nwmx),br(nwmx,1)
      real      dt

      complex wrk(nwmx),rnrmt(nchmx,ndmx,nt),cwrk(nchmx,nsmax)
     &   ,cout(nchmx,nsmax)
        
      integer ix1(nchmx1,nmax),
     &   
     &   next(ndmx),ifirstd(ndmx),nacc(ndmx),
     &    idp(0:ndmx,2,nbadmx),iwrk(nchmx,nsmax),
     &    ids(0:ndmx,2,nbadmx),ndp(0:ndmx),ifirst(ndmx),nspec(ndmx)
     4   ,nds(0:ndmx),iftype(nfilmax,nchmx),
     5    nfil(nchmx),iset(4,nsmax),ifdir(3,nsmax) 
     6    ,idl(ndmx),ixs(nsmax)
      integer iset_offset(ndmx)

        logical lstart(ndmx),lstartd(ndmx),lgood,lg,lclkd
     &  ,lfd(nchmx,ndmx),lwinst(ndmx),lend

      character*1 ans  ! c0,
      character*50 arg
      character*6 chid(nchmx)
      character*80 afparam(nfilmax,nchmx)

C   SCRATCH FILE 02.01.98
      real x_scr(nsmax*nfmax,nchmx,2)
      integer iuse_scr(nsmax*nfmax)

       print*,'********************************************************'
       print*,'*                                                      *'
       print*,'*     DNFF     Version: ',cvers,'   Date: ',cdate,
     &'        *'
       print*,'*                                                      *'
       print*,'********************************************************'
       print*,' '
       bytes = NBYTES
       l_asc = (BINASC .eq. 'ASC')
       l_clk_hd = .false.

ceg    ioff = 0
ccc      (ioff is offset to add to sample numbers returned by reading
ccc       routine;  this is set in bdrcsu)

c     parse command line options
      nargs = iargc()
      do k = 1, nargs
         call getarg(k,arg)
         if (arg(1:2) .eq. '-b') then
            read (arg(3:3),*) bytes
         elseif(arg(1:2) .eq. '-a') then
ccc         ascii data file, read clock reset info from .clk file
            l_asc = .true.
         elseif(arg(1:2) .eq. '-A') then
ccc         ascii data file, read clock reset info from header of data file
            l_asc = .true.
            l_clk_hd = .true.
         endif
      enddo

      if(l_asc) then
         write(6,*) '  Reading from ASCII data files'
         write(6,*)
         if(l_clk_hd) then
            write(6,*) '  Assuming clock reset information ',
     &                  'in data file headers'
            write(6,*)
            write(6,*)
         else
            write(6,*) '  Reading clock reset information from ',
     &                  '.clk file'
            write(6,*)
            write(6,*)
         endif
      else 
         write(6,*) '  Reading from binary data files'
         write(6,
     &     '(/,''   Using integer*'',i1,'' as binary input '',/)')
     &     bytes
      endif

ccc   return here to start processing another data file
1000  continue

c*********************************************************************

c    compute l1 power spectrum for pre-scaling output
c          blank array pspecl1
      do i = 1,ndmx
         nspec(i) = 0
         do k = 1,nchmx
            do j = 1,nwmx
               pspecl1(k,j,i) = 0.0
            enddo
         enddo
      enddo

c*****  get station name from standard input, input files paths from
c       local file paths; open input file, make sp br etc path/file names
      call cininit(msval,    !cfsp,inunit,
     &    rmsval,nch)        !,nch

      nchp1 = nch+1
      call outinit(nch,lpack)

c  read system parameter file:  number of channels,
c  electrode line lengths, filter parameters, conversion
c  factors, etc. from file spsta###x
      call getsp(nch0,sampr,sc,nfil,iftype,afparam,
     &   decl,stcor,orient,chid,cda,cdb,lclkd)

      if(nch0.ne.nch) then
         print*,'error: nch in sp file does not agree with nch in',
     &     ' data file'
         stop
      endif
      if(cfsp(4:11).eq.'standard') then
         print*,'Enter declination, station coordinates'
         read(5,*) decl,stcor(1),stcor(2)
      endif

c  set up decimation; decset reads parameters which control decimation
c  from a file; see documentation in decset for further info; 
      call decset(nwmx,lfd,idoff,ierr)
      if(ierr.eq.1) then
         print*,'parameters read from cf_decset are inconsistent',
     &   ' with parameters declared in main program'
           stop
      endif

c  compute sampling rate for each decimation level     
        nfusemx = 0
        do 1 i = 1,nd
        idl(i) = i + idoff
        samprate(i) = sampr
        nfusemx = max(nfusemx,nfuse(i))
           do 1 j = 1,i                        
1          samprate(i) = samprate(i)*idec(j)
c  set up decimation filter corrections, conversion of counts
c  to physical units (this routine makes a table of transfer
c  functions to correct for all of these)
      call fcorsu(samprate,nfil,afparam, !fc,nfc,nfcmx,nwin,nd,nch,
     1  iftype,sc,rnrmt,cda,lfd)
c    read in bad records file (gives record ranges which should not 
c    be processed (idp) or stacked (ids); see routine for input format
c    file also gives offset to adjust record numbers with
        call bdrcsu(nd,nbadmx,ndp,idp,nds,ids,ioff)
ccccEGBERT
ccc    Change call to mk_offset, and replace mk_offst.f with new version
ccc     declare dt (output by mk_offset, and used by phs_shft) as real
           call mk_offset(time_int,iset_offset,sampfreq,dt)
ccccEGBERT
c      write(*,*) 'TIME_INTERVAL',time_int

c    set up fft tables
      do 5 id = 1,nd
         call cffti(nwin(id),wsave(1,id))
5        continue

c  initialize various arrays
      do 6 i = 1,nd
6        lwinst(i) = .true.
      nuse = 0
      nmsmx = nwmx

c    get first block of data and set up pointers at start of processing
c    also, if too large a data gap occurs return to stmt. 90
c    and reset pointers


      lfirst = .true.
      if(l_asc) then
         call rdasc(npts,ix1,ierr,nch)
      else
         call rdblk(npts,ix1,ierr,nch)
      endif
c      write(0,*) 'npts',npts
      if(ierr.eq.-1) go to 300
      lend = .false.
90    irec1 = ix1(1,1)
      call pterst(irec1,ifirst,ifirstd,next,lstart,lstartd)

c>>>>>>>>>>>>>>>main loop starts here >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ifrec = 0
100     continue
        
c     call decimation routine
        call dcimte(ix1,npts,xx,nchp1,nmax,ifirst,ifirstd,next,
     &     lstart,lstartd,nacc,rmsval)
                           
      do 200 id = 1,nd
      if(nfuse(id).gt.0) then
c   compute number of sets accumulated at decimation level id
c    (will ususally be 1 or 0)
      nsets = (nacc(id) - olap(id) - 2*npwmx(id) - missmx(1,id)/2)
     &                           / ( nwin(id) - olap(id) )
         do 190 i = 1,nsets                                  
c          make set i at decimation level id; return in array x
           call mkset(xx,nchp1,nmax,ifirst,x,nwmx,ist,id,irec,
     &         lrec,rmsval,lgood,roff)

            ist = ist + iset_offset(id)
c   need to decide if set contains data segmenta which should not be
c   processed (declared in bad record file brsta###x)
           call badrec(nd,ndp,idp,id,irec,lrec,lg)
           if (lgood.and.lg) then
c       set is "good" (i.e. not too many missing data points; not
c       among the bad record ranges); proceed with processing

c   first difference,demean, window, and fourier transform data 
             n = nwin(id)
             nmen=n+npwmx(id)-1
             call demean(x,nch,nmen)
             call frstdif(x,nch,n,lfd,id,areg,npw(1,id),rr,ar,br,bwo)
             call demean(x,nch,n)
             call ps1win(x,nch,n,lwinst(id),w(1,id))
             call fftnp2(n,nch,x,wrk,wsave(1,id))

c   correct for filters (decimation, first differencing) and normalize
c   fourier coefficients to have units of nt/(sqrt(hz))
c   [ uses table of transfer function coefficients stored in complex
c     array rnrmt ]
            nfmid = nfuse(id)
            call filtcor(x,nch,nfmid,id,rnrmt,nd,areg,npw(1,id),n,lfd)
            nfreq = nfmid

ccccEGBERT
ccc ADD THIS CALL ... and include new subroutine phs_shft.f
cc      correct for offsets in actual sampling times
           call phs_shft(dt,x,nfreq,nch,samprate(id),nwin(id))
ccccEGBERT

c        if necessary correct for clock drifts
              if(lclkd) then          
                 tset = ((irec+lrec)/2. - iclk0)*samprate(1)
                 call cldrft(x,nfreq,nch,samprate(id)
     &              ,nwin(id),cdb,tset)
              endif

c     if necessary correct for offset of begining of set (caused by
c        non-integer set overlap)
              if(abs(roff) .gt. 1.0e-5) then
                 roff = - roff * samprate(1)
                 call cldrft(x,nfreq,nch,samprate(id)
     &              ,nwin(id),1.0,roff)
              endif

c       check to see if set contains segments which should not be stacked
c       if so indicate by making set number negative
              call badrec(nd,nds,ids,id,irec,lrec,lg)
              if( .not. lg) ist = - abs(ist)

c         output fourier coefficients to scratch file f_sta###x
c              call freqout(ioscr,nch,nfreq,x,iuse,ifrec)
c              print*,'CALL FREQOUT: NFUSEMX', nfusemx
              call freqout(nch,nfreq,x,ifrec,nsmax,nfusemx,  ! SCRATCH FILE 02.01.98
     &  x_scr,iuse_scr)

c         add abs(FCs) to array to compute scalling factor for output
              call l1spec(nch,x,nfreq,pspecl1(1,1,id),nchmx)
              nspec(id) = nspec(id) + 1

c         update array ifdir (which will be used to sort out fourier coefficients
c         at end of program)
              nuse = nuse + 1
              ifdir(1,nuse) = id + idoff
              ifdir(2,nuse) = nfreq
              ifdir(3,nuse) = ist

c     check for overflow of nsmax (maximum # of sets
              if (nuse.eq.nsmax) then
                 write(*,232) nsmax
 232   format(/,'!!!  Number of processed sets reached',/,
     &          '     maximum number of allowed sets !!!',/,
     &          '     The program will finish processing now',/,
     &          '     To process all data change NSMAX =',i6,/,
     &          '     in file PARAMS1.INC !!!')
                 goto 300
              endif

              print*, 'good',irec,lrec,id,ist
c           else
c            print*,'bad',irec,lrec,id,ist
           endif
190     continue
        endif
200     continue
c    get more data and go back to top of main loop
      if(l_asc) then
ccc      slight complication here to get exactly the same results
ccc      with both ascii and binary.  With ASCII file EOF is detected
ccc      in rdasc after reading some data.  Set lend upon finding EOF,
ccc      use the data from the final read to make more sets, then
ccc      terminate.  With binary files EOF is not detected until
ccc      after the last data block has been read so this comlication
ccc      is not necessary
         if(lend) then
            goto 300
         else
            call rdasc(npts,ix1,ierr,nch)
            lend = (ierr .eq. -1)
            ierr = 0
         endif
      else
         call rdblk(npts,ix1,ierr,nch)
      endif
                    
c     ierr = 0 for normal read
        if(ierr.eq.0) then
           lfirst = .false.
           go to 100
        endif

c      ierr = -3 when there is a gap in the record numbers that is so
c large that it is most reasonable to start the decimation filtering over
        if(ierr.eq. -3) then
           lfirst = .true.
           go to 90
        endif

c     ierr = -2 means the record numbers jump backwards; the program doesn't
c    know what to do about this and stops
        if(ierr.eq.-2) then
            print*,'error: excess records/ incorrect time mark'
            print*,'stopping here'
        endif

300     continue

c<<<<<<<<<<<<<<<<<<< end of main loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        print*
      print*,'number of sets = ',nuse
        print*
c***********************
c       reorder f_ file by frequency
c*************************

c       close connection to data file
      close(in_unit)
c       reorder and output file of fourier coefficients
c      l1 power spectrum
      do 490 i = 1,nd
      do 490 j = 1,nch
      do 490 k = 1,nfuse(i)
         if(lpack) then
            pspecl1(j,k,i) = pspecl1(j,k,i)/nspec(i)
         else
            pspecl1(j,k,i) = 1000.
         endif
490      continue

ccc   eliminate decimation levels with nfuse(id) = 0
      iid = 0
      do id = 1,nd
         if(nfuse(id) .gt. 0) then
            iid = iid + 1
            idl(iid) = idl(id)
            nwin(iid) = nwin(id)
            samprate(iid) = samprate(id)
         endif
      enddo

c      header for FC file
      if(lpack) then
         irlo = 4*(nch+1)
      else
         irlo = 4*(2*nch+1)
      end if

      call wfhead(nch,iid,nfusemx,nwin,samprate,idl
     1   ,chid,orient,decl,stcor,irlo)

c       read in scratch file, reorder pack if requested, output        
C     ME        REMOVE SCRATCH FILE
c      call mkseq(ioscr,iouout,ifdir,iset,nuse,nch,nd,
c     1 cwrk,cout,nfusemx,idoff,pspecl1,nchmx,nwmx,lpack,iwrk,ixs)
      call mkseq(ifdir,iset,nuse,nch,nd,
     1 cwrk,cout,nfusemx,idoff,pspecl1,nchmx,nwmx,lpack,iwrk,ixs,
     2 nsmax, x_scr, iuse_scr)
C     ME        REMOVE SCRATCH FILE        
c       close files; delete f_file on closing
      close(out_unit)
C ME    REMOVE SCRATCH FILE
c      close(ioscr,status = 'delete')
C ME    REMOVE SCRATCH FILE
c500   continue

      write(0,'(''FT another data file?  '',$)')
      read(5,'(a1)') ans
      if(ans.eq.'y') go to 1000
      stop
      end
