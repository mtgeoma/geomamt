        program tranrr
        character cvers*3, cdate*10
        parameter (cvers = '2.6', cdate = '03/10/1998')
c       MT transfer function program; will handle remote reference
c         etc.

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c     array size parameters are set in iosize.inc
      include 'iosize.inc'

c     parameter tol gives tolerance for testing convergence of
c     robust algorithm - this is the fractional change in the tf
c     between iterations
c       COHW controls relative width of coherence presort band
c       ICOHBM is minimum width of coherece presort band
      parameter( tol = .002,cohw=1.5,icohbm=10)
c*********************************************
c     03/10/98        ME
c
c     command line options added
c     -x output sdm in file
c     -s read new channel ordering from grp-file
c     -p print out first/last set numbers of each FC file
c     -S process only sets within specified windows
c
c     output to console cleaned 
c     
c*********************************************
c===================================================================

      integer iband(2,nbmax),idl(nbmax),ldf(nbmax),
     &  nch(nstamx),ntape(nstamx),nd(nstamx),npts(ndmax,nstamx),
     &    ibcoh(2),ibt(2),iwrk(nstamx,ntfmax),
     &   ixs(nfull),ixf(nfull)
     
      real dr(ndmax,nstamx),stcor(2,nstamx),
     &        orient(2,nchmx,nstamx),decl(nstamx),period(nbmax),
     & w(ntfmax),rot(20),rdf(ntmx,nbmax),coh(3,nsetmx),cohp(2)
 
      logical YoN,lpcoh,lrbst,lbuse(nbmax),lprint,lrref,leref,lcoh,
     & lx(nstamx,nlx)
        
      complex xx(ntmx,ncbmx),xo(ntmx,ntfmax),   !zt(nsmx,nbmax),
     &   z(nsmx,nbmax),zzdum(nsmx,nbmax) !  ,zz1(nsmx,nbmax),zz2(nsmx,nbmax)
        
      character*40 cbandf,cfpcoh,stnames(nstamx,ntpmax)
      character*80 outname
      character*20 cfile
      character*6 chid(ntmx)
      character*1 ans
      character*80 chead
cnew      character*60 ctitle(20)
      character*30 cdirout,cdirin
      character*80 arg                 ! command line arg

cnew variables for reordering
      logical         l_GRP, l_SDM, l_PRISET, l_SETS
      character*80    cfile_GRP
      character*6     tempc(ntmx)
      integer         iin(2), iout(ntmx), irref(2) ,nchout
      real            tempr(2,nchmx,nstamx)
      integer         unit_GRP
      PARAMETER (unit_GRP = 17)

c ME
      integer         isuse(2,ndmax,10), nuset, dec_fac
c ME

      common /options/lrbst,lrref,leref,lcoh,coht,cohp,cohm,num
     & ,nrot,rot,cdirout,cdirin,chead,cfpcoh,lpcoh

      YoN(ans) = (ans.eq.'Y').or.(ans.eq.'y')

c ME

      print*,' '
      print*,'*********************************************************'
      print*,'*        TRANMT   Version: ',cvers,' Date: ',cdate,'    *'
      print*,'*********************************************************'
      print*,' '
      print*,'    OUTPUT: Z-file format. Extensions zss (single site),'
      print*,'    zrr (remote refernce) are appended to the name'
      print*,'    specified in tranmt.cfg.'
      print*,'    '
      print*,'    Command line options: '
      print*,'      -s<group.cfg> for non-default grouping of channels'
      print*,'      -x output spectral density matrix (.sdm) '
      print*,'      -p print out # of first/last set per site and stop '
      print*,'      -S<sets> process only specified sets  '
      print*,'             (see documentation for details).'
      print*,' '
c ME

c     lprint is true to print out debugging statements
      lprint = .false.
      l_GRP = .false.
      l_SDM = .false.
      l_PRISET = .false.
      l_SETS = .false.
c set this to make sure all sets are processed if no -S option is specified
            maxint = 2147483647
            nuset = 1
            isuse(1,1,1) = - maxint
            isuse(2,1,1) = maxint

c*******************************************************************
        nargs = iargc()
        do iargs = 1, nargs
           call getarg(iargs, arg)
           if (arg(1:2) .eq. '-s' ) then
              l_GRP = .true.
              cfile_grp = arg(3:80)
           elseif (arg(1:2) .eq.'-x') then
c     write SDM in intermediate GEOTOOLS format to file
              l_SDM = .true.
           elseif (arg(1:2) .eq.'-p') then
c     only print out numbers of first and last set of each site
c     to set up set limits to be processed
              l_PRISET = .true.
           elseif (arg(1:2) .eq.'-S') then
              l_SETS = .true.
              read(arg(3:80),*) 
     &             nuset,(isuse(1,1,ka),isuse(2,1,ka), ka = 1,nuset)
c              print*, nuset, (isuse(1,1,ka),isuse(2,1,ka), ka = 1,nuset)
           endif
        enddo


      do i = 1,nbmax
         ldf(i)=0
         lbuse(i) = .true.
      enddo




      print*,'enter control file name'
      read(5,'(a20)') cfile
      print*
      open(unit=99,file = cfile)

c>>>>>>>>> return here to process another data set
9999  continue

c.....get processing info, open files for io, set up bands, etc.
c     ME        isuse, and nuset is returned from  SETUP
      call setup(nsta,ntape,stnames,outname,nd,nch,cbandf,stcor,decl,
     &orient,chid,dr,npts,nchstk,ncht)  !,isuse,nuset


c.... set up for frequency band averaging
      call bset(nbt,iband,idl,period,ndmax,nbmax,ierr,
     &   dr,npts,cbandf)   !,isuse,nuset

c     ME        Compute the sets to use for processing
c      if (nuset.ge.1) then
c         nduset = 0
c         do ista = 1, nsta
c            nduset = max(nduset,nd(ista))
            nduset = idl(nbt)
c         enddo
         do i = 1, nuset
            do k = 2, nduset
               dec_fac = nint(dr(k,1)/dr(k-1,1))
            isuse(1,k,i) = int(isuse(1,k-1,i)/dec_fac)
            isuse(2,k,i) = int(isuse(2,k-1,i)/dec_fac)
            if (mod(isuse(1,k-1,i),dec_fac).gt.0) 
     &           isuse(1,k,i)=isuse(1,k,i)+1 
            if (mod(isuse(2,k-1,i),dec_fac).gt.0) 
     &           isuse(2,k,i)=isuse(2,k,i)+1 
            enddo
         enddo 
c         if (isuse(1,1,1).ne.-maxint .and. isuse(2,1,1).ne.maxint) then
         if (l_SETS) then
            write (6,894)
            write (6,890)
            write (6,891) (i,i, i=1,nuset)
            write (6,892)
            do k = 1, nduset
               write (6,893) k, (isuse(1,k,i),isuse(2,k,i), i = 1,nuset)
            enddo 
            write (6,894)
         end if

 890     format ('              SET-WINDOWS FOR PROCESSING ')
 891     format ('SET-WINDOW:',1x,5(i2,'/1',2x,i2,'/2',2x))
 892     format ('DEZ-LEVEL')
 893     format (i5,5x,20i6)
 894     format (/,'***********************************************',/)
c895     format (/,'                 PROCESSING ALL SETS  ',/)

c      endif
c     ME        end of computing isuse

 
c.....make directory for fourier coefficient files; open files
      call mkfdir(nsta,ntape)

      if(.not.lfop) then
         call fop(nsta,nch,ntape)
      end if

cnew if reordering of channels is requested, read new channel setup
cnew from cfile_GRP
      if (l_GRP) then
         print*
         print*,'Reading new channel grouping from ',
     &        cfile_grp(1:iclong(cfile_grp,80))
         open (unit=unit_GRP, file=cfile_GRP, status='old', err=3001)
         read (unit_GRP,*) iin
         read (unit_GRP,*) nchout
         read (unit_GRP,*) (iout(k), k = 1, nchout)
         nchstk = nchout + 2

         write(*,'(/,'' New predicting channels:'',2i4)') iin
         write(*,'(27x,a6,2x,a6)') (chid(iin(l)),l = 1,2)
         write(*,'(/,'' New predicted channels: '',18i4)') 
     &        (iout(k),k=1,nchout)
         write(*,'(27x,18(a6,2x))') (chid(iout(l)),l = 1,nchout)

         if (lrref) then
            read (unit_GRP,*) irref
            nchstk = nchout + 4
         write(*,'(/,'' New reference channels:  '',2i4)') irref
         write(*,'(27x,a6,2x,a6,/)') (chid(irref(l)),l = 1,2)

         endif
         close (unit_GRP)
      endif
*******************************************************************
c**************** LOOP EXECUTED ONCE FOR EACH FREQUENCY BAND *******
c*******************************************************************
      write(*,*) 
      write(*,*) 
     & 'IBAND D-Lev   Period   IBAND(1) IBAND(2) INITIAL FINAL ITER'
      do 100 ib=1,nbt
         if(lprint) then
            print*,'frequencies',iband(1,ib),'  to   '
     1                             ,iband(2,ib),'level ',idl(ib) 
            print*,'period',period(ib)
         end if
         id=idl(ib)
c........find limits of coherence pre-sort band
         if(lcoh) then
            call cohband(iband(1,ib),ibcoh,nfreqmax,nfb,ncb,cohw,
     &                                                     icohbm)
            print*,'coherence band: ',ibcoh
         else
            ibcoh(1) = iband(1,ib)
            ibcoh(2) = iband(2,ib)
            nfb = iband(2,ib)-iband(1,ib)+1
            ncb = nfb
         end if

c........make array containing all data for coherence sort band
c            (includes data for narrower estimation band)

c       ME      nuset added to pass to mkrec
         call mkrec(nd,id,ibcoh,nsta,ntape,nch,ncht,isuse,nuset,nfreq,
     &         xx,lx,ixs,ixf,xo,iwrk,l_PRISET)

cnew if reordering is requested, reorder channels now
         if (l_GRP) then
            call reorder(xx,ncht,nfreq,iin,nchout,iout,lrref,irref)
         endif
         
c..........coherence presort
         ns1 = nfreq/(ibcoh(2)-ibcoh(1)+1)
         ncb = nfreq/ns1
         if(lcoh) then
            ibt(1) = iband(1,ib) - ibcoh(1) + 1
            ibt(2) = iband(2,ib) - ibcoh(1) + 1
            ibt(2) = min(ibt(2),ncb)
            call cohsrt(xx,ncht,lrref,ncb,ns1,ibt,cohm,num,coht,
     &       cohp,nu,coh,lpcoh,cfpcoh,xo,period(ib),w,ntfmax)

            if(lprint) then             
               print*, 'data points in coh band',nfreq
               print*, 'no. of data points used',nu
            end if
         end if

         if(lcoh) then
            call rxspclev(xo,nu,nchstk,z(1,ib),ldf(ib),rdf(1,ib),tol
     &      ,lrbst,w,ncht,lrref,iiter)
         else
            call rxspclev(xx,nfreq,nchstk,z(1,ib),ldf(ib),rdf(1,ib),tol
     &      ,lrbst,w,ncht,lrref,iiter)
         end if
         write(*,900) ib,idl(ib),period(ib),
     &        iband(1,ib), iband(2,ib),nfreq,ldf(ib),iiter
 900     format(1x,i4,3x,i3,1x,1pe12.6,1x,i6,3x,i6,3x,i5,2x,i5,1x,i4)
100     continue

      close(2)
      close(21)
      if(lrref) close(22)
      write(*,*) 
     & 'IBAND D-Lev   Period   IBAND(1) IBAND(2) INITIAL FINAL ITER'

c************END OF FREQUENCY BAND STACKING LOOP *********************
ccc   output transfer functions with error covariance info
      ns = (nchstk*(nchstk+1))/2
ccc   reordering of orientations and channel id's to reflect
ccc   user supplied channel order
      if(l_GRP) then
        tempc(1)=chid(iin(1))
        tempc(2)=chid(iin(2))
        do j=1,2
          tempr(j,1,1)=orient(j,iin(1),1)
          tempr(j,2,1)=orient(j,iin(2),1)
        enddo
        do i=3,nchstk-2
          tempc(i)=chid(iout(i-2))
          do j=1,2
            tempr(j,i,1)=orient(j,iout(i-2),1)
          enddo
        enddo
        tempc(i+1)=chid(irref(1))
        tempc(i+2)=chid(irref(2))
        do j=1,2
          tempr(j,i+1,1)=orient(j,irref(1),1)
          tempr(j,i+2,1)=orient(j,irref(2),1)
        enddo
        do j=1,nchstk
           chid(j)=tempc(j)
           do i=1,2
             orient(i,j,1)=tempr(i,j,1)
           enddo
        enddo 
      endif
c      write(*,'(10a6)') (chid(k),k=1,ncht)
       

      call wrt_z(z,nbt,period,rdf,ldf,outname,orient,
     +   nchstk,ns,lrref,nsmx,zzdum,iband,cdirout,idl,
     +   stcor,decl,chid,dr,ndmax,chead,ntmx,l_SDM,bw)

c      print*,'do another station?'
      read(99,'(a1)') ans
      if( YoN(ans)) go to 9999
C200   continue
      close(99)
      stop
 3001 print*,'Specified GROUP file not found!'
      STOP

      end
