c_______________________________________________________________
c
c ME   02/26/97 added isuse 3-d array of sets to use for processing to be passed
c ME		to SETUP. See comments marked ME below

      subroutine setup(nsta,ntape,stnames,outname,nd,nch,cbandf,stcor,
     &  decl,orient,chid,dr,npts,nchstk,ncht) !,isuse,nuset

      include 'iosize.inc' 
      include '../include/four_byte.inc'
      character*6 chid(*)
      character*1 ans
      integer nch(*),ntape(*),npts(ndmax,*),nd(*)
      character*40 stnames(nstamx,ntpmax),fname,cfpcoh
      character*30 cdirin,cdirout !,cdircf
      character*80 chead,cfile,outname,cbandf

      logical YoN,lrref,leref,lpcoh,lcoh,lrbst
      real dr(ndmax,*),orient(2,*),stcor(2,*),decl(*),rot(20),cohp(2)

      common /options/lrbst,lrref,leref,lcoh,coht,cohp,cohm,num,
     &nrot,rot,cdirout,cdirin,chead,cfpcoh,lpcoh

      YoN(ans) = (ans.eq.'Y').or.(ans.eq.'y')

c      open unit 1 before calling

c..........read in output identifier
      read(99,'(a80)') outname
c..........read in options file; open and read
      read(99,'(a80)')  cfile
      open(unit=2,file=cfile)
c......... read in header describing options
      read(2,'(a80)') chead
c        directories: input files
c        directories: output files
      read(2,'(a30)') cdirin
      read(2,'(a30)') cdirout
c          bandset up file
      read(2,'(a80)') cbandf

c      print*, 'robust ? (y/n)'
      read(2,'(a1)') ans
      lrbst = YoN(ans)

c      print*, 'remote reference ? (y/n)'
      read(2,'(a1)') ans
      lrref = YoN(ans)

c      print*, 'electric field reference ? (y/n)'
      read(2,'(a1)') ans
      leref = YoN(ans)

c      print*,'output coherence (y/n)?'
      read(2,'(a1)') ans
      lpcoh = YoN(ans)
      if(lpcoh) then
c         print*,'enter file for coherence output'
         read(2,'(a20)')  cfpcoh
      end if

c      print*,'enter coherence sorting parameters'
c         coht - is target coherence; ideally would like
c          everything to achieve this coherence
c         cohp - are used to determine target number of points
c          via : nutarget = cohp(1)*nu**cohp(2) where nu is number
c          of points available (e.g.: cohp(1) = 3., cohp(2) = .5)
c          If possible would like to get nutarget data points to
c          use; if there are not nutarget points in sets with
c          coherence above coht, accept lower coherence
c          sets, until ntarget points are available, or until:
c         cohm - the minimum acceptable coherence is reached
c          never accept sets with lower coherence unless:
c         num - the absolute minimum number of points can't
c          be found which are in sets achieving a coherence of cohm.
c
      read(2,*) coht,cohp,cohm,num
      lcoh = (coht .gt. .00)

c..........read in rotation angles
c      read(2,*) nrot
c      read(2,*) (rot(i),i=1,nrot)

c*****done reading options file
      close(2)

c     Data file info
c     number of stations 
      read(99,*) nsta

c     for each station, number of tapes, channels
      do i = 1,nsta
         read(99,*) ntape(i),nch(i)
ccc      for each station-tape, FC file name
         do j = 1,ntape(i)
            read(99,'(a40)') stnames(i,j)
         enddo
      enddo


      ncht = 0
      do i = 1,nsta
         ncht = ncht+nch(i)
      enddo

      if(lrref) then
         nchstk=ncht-nch(nsta)+2
      else
         nchstk=ncht
      end if


ccc   open data files
      m = iclong(cdirin,30)
      ii=20
      ista = 1
      do i=1,nsta
         if(lpack) then
            iorecl(i) = 4*(nch(i)+1)
         else
            iorecl(i) = 4*(2*nch(i)+1)
         end if
         if(l_4byte) iorecl(i) = iorecl(i)/4

         print*,'station # ',i
         do j=1,ntape(i)
            ii=ii+1
            inunit(i,j)=ii
            fname=stnames(i,j)
            if(m.gt.0) then
               cfilein(i,j) = cdirin(1:m)//'/'//fname
            else
               cfilein(i,j) = fname
            end if
c            write(*,*) cfilein(i,j)
            open(unit=ii,file=cfilein(i,j),
     &         form='unformatted',access='direct',recl=iorecl(i))
            write(*,*)'unit ',inunit(i,j),'  file  ',
     &        cfilein(i,j)(1:iclong(cfilein(i,j),80)),
     &           '  recl=',iorecl(i)
         enddo      ! j=1,ntape(i)

ccc      read header for first tape for each station ... 
ccc       get channel IDs, orientations etc.
         call rfhead(inunit(i,1),nchi,ndmax,nd(i),nfmaxi
     &    ,npts(1,i),dr(1,i),chid(ista),orient(1,ista),decl(i),
     &          stcor(1,i),iorecl(i))
         ista = ista + nch(i)

         if(nchi.ne.nch(i)) then
            print*,'error: nch in set up file disagrees with nch in'
     &                      ,' file header'
            print*,'file = ',fname,' nch in set up file = ' 
     &            ,nch(i), '  nch in header = ',nchi
            return
         end if
         if(nd(i).gt.ndmax) then
            print*,'error: nd(i) .gt. ndmax','  i = ',i,' nd= ',nd(i),
     &         'ndmax =',ndmax
            return
         end if
         if(nfmaxi.gt.nfreqmax) then
            print*,'error: nfreqmax in file header exceeds max set',
     &       ' in program'
            print*,'file = ', fname,'nfmaxi =',nfmaxi,
     &               '  nfreqmx = ',nfreqmax
         end if
      enddo     ! i = 1,nsta

ccc      count number of E channels (excluding remote site ...)
         if(lrref)  then
            nsta1 = nsta - 1
         else
            nsta1 = nsta
         endif
         nech = 0 
         i1 = 0
         do ista = 1,nsta1
            do ich = 1,nch(ista)
               i1 = i1 + 1
               if(chid(i1)(1:1) .eq. 'E') nech = nech+1
            enddo
         enddo

      return
      end
c__________________________________________________________________
c
      function iclong(croot,n)
      character*1 croot(n)
      do 10 i = n,1,-1
      if(ichar(croot(i)).ne.32) then
         iclong = i
         return
      end if
10    continue
      iclong = 0
      return
      end

