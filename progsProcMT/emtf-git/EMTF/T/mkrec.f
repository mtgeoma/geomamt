c___________________________________________________________________
      subroutine mkrec(nd,id,iband,nsta,ntape,nch,ncht,isuse,nuset,
     &  nfreq,xx,lx,ixs,ixf,cwrk,iwrk,l_PRISET)
 
c  subroutine makes complex FC array records (stored in array xx) for ntape
c   tapes of data, nsta stations; record includes frequencies iband(1)
c   through iband(2) at decimation level id; requires the directory array
c   irecd (in commonblock ioblock) be constructed by a call to
c      routine mkfdir bdfore first call to mkrec;
c    other inputs: ncht = total channels of data
c     isuse: Only sets at decimation level id with numbers between
c         isuse(1,id) and isuse(2,id) will be output


c      logical variable llx is a data item which controls how the routine
c       handles sets for which some stations do not have data      
c    IF llx = .true.  the routine returns a logical array lx,
c       the elements of which are true or false depending on whether
c       good data is present for a given station at a given frequency;
c    In this case storage for the logical array lx of total dimension
c     nsta*nfreq must be provided for in the calling program;
c    if llx = .false. lx only has to be of dimension nsta
c      in this case only sets for which all data are present and good
c      will be returned; standard usage for routine TF processing

c     if lfull = .true. return set number (in ixs) and frequency band
c       number (in ixf); otherwise just FCs

c     If lfop is .true. open and close files for each call to readfvg
c      use this way when number of files is large

      include 'iosize.inc'

      logical lx(nsta,*),lset(100),luse,lusethis, l_PRISET
c     ME        isuse array changed
      integer ntape(*),
     1  iset(100),lstset(100),iband(2 ),ncht0(100),iwrk(nsta,*),
     2  isuse(2,ndmax,10),nch(nsta),ixs(*),ixf(*)
      complex cwrk(ncht,*),xx(ncht,*)

      data maxint/2147483647/

c******* initialize
      ifreql = 1
      ifreq = 1
      ncht0(1) = 0
      do 1 i = 1,nsta-1
         ncht0(i+1) = ncht0(i) + nch(i)
1        continue

c*>>outer loop: collect and sort data for each frequency in estimation band
      do 100 ifb=iband(1),iband(2)
         nmax=0
c********* read loop for nsta stations
          do 10 ista = 1 , nsta
             i0=1
             do 5 itape = 1,ntape(ista)
                if (lfop) then
c                 need to open file before reading
                   open(unit=inunit(ista,itape),
     &             file=cfilein(ista,itape),form='unformatted',
     &             recl=iorecl(ista),access='direct')
                end if
                call readfvg(id,ifb,ista,itape,nd,nsta,ntape,
     &            nch(ista),cwrk(ncht0(ista)+1,i0),
     &            iwrk(1,i0),nf,ncht)
                if(lfop) close(inunit(ista,itape))
                i0 = i0 + nf

c     ME if print out of set nnumbers is requested (command line option
c     ME -p) print out iwrk(ista,1),iwrk(ista,1)+nf 
                if (l_PRISET) then
                   write(*,'(''Station '',i3,'' FILE '',i3,
     &                  '' FIRST/LAST SET '', 2i6)')
     &                  ista,itape,iwrk(ista,1),iwrk(ista,1)+nf 
                endif
5               continue
             nf = i0 - 1
c               nmax is maximum number of sets over all stations
             if (nf .gt. nmax) nmax = nf
c                 mark the end of set number array 
             iwrk(ista,nf+1)=maxint
c                initialize iset, lstset
             do 7 k = 1,nf
c     ME        isuse array changed
                if(abs(iwrk(ista,k)).ge.isuse(1,id,1)) go to 8
7               continue
             k = nf+ 1
8            continue
             iset(ista)=abs(iwrk(ista,k))
             lset(ista) = (iwrk(ista,k).ge.0).or.ljunk(ista)
             lstset(ista)=k
10           continue

c     ME if only print out of set numbers is requested, stop here
             if (l_PRISET) stop

c*********  data is now read in for frequency ifb, and stored in array cwrk
c           Now: organize into records with all set numers the same for
c               all stations
 
c>>>>>>>>>>>>>>>>>>>>>> START OF SORTING LOOP >>>>>>>>>>>>>>>>>>>>>>>>>>>
15         continue
c            find smallest set number
           nxtset=iset(1)
           do 20 ista=2,nsta
20            if(nxtset.gt.iset(ista)) nxtset=iset(ista)
 
c******** if nxtset .eq. maxint  then done with frequency ifb
c     ME        isuse array changed
           if((nxtset.eq.maxint).or.(nxtset.gt.isuse(2,id,nuset)))
     &             go to 100
 
c     ME +++++++++++++++++++++
c     if nxtset is not within the specified boundaries of sets to be processed
           lusethis = .false.
           do iuset = 1,nuset
              if ((nxtset.ge.isuse(1,id,iuset)).and.
     &                  (nxtset.le.isuse(2,id,iuset))) then
                  lusethis = .true.
              end if
           enddo 

           if (.not.lusethis) then
              do ista = 1, nsta
                 lstset(ista) = lstset(ista) + 1
                 iset(ista) = abs(iwrk(ista,lstset(ista)))
                 lset(ista) = (iwrk(ista,lstset(ista)).gt.0)
     &                .or.ljunk(ista)
              enddo 
              goto 15
              
           endif
c     ME ++++++++++++++++++++++

c****** else put together record for set number nxtset
c       xx contains complex FCs
c      AND (optionally)
c       ixs contains set numbers (if lfull = .true.)
c       ixf contains frequency number (ifb) (if lfull = .true.)
c       lx is logical array, true if data is present and OK for station/set
c                            false if not (if llx = .true.)
          if(lfull) then
             ixs(ifreq)=nxtset
             ixf(ifreq)=ifb
          end if
          do 30 ista = 1,nsta
             if (iset(ista).eq.nxtset) then
                 do 25 i=1,nch(ista)
25                  xx(ncht0(ista)+i,ifreq)
     &                        =cwrk(ncht0(ista)+i,lstset(ista))
                 lx(ista,ifreql)= lset(ista)
                 lstset(ista)=lstset(ista)+1
                 iset(ista)=abs(iwrk(ista,lstset(ista)))
                 lset(ista)=(iwrk(ista,lstset(ista)).gt.0)
     &                                       .or.ljunk(ista)
              else
                 lx(ista,ifreql) = .false.
              end if
30            continue
           if(llx) then
              ifreql = ifreql + 1
              ifreq=ifreq+1
           else
              luse = lx(1,1)
              do 35 ista = 2,nsta
35               luse = luse.and.lx(ista,1) 
              if(luse) ifreq = ifreq+1
           end if

c            check to see if last frequency has been used
           if(ifreq.eq.ncbmx) go to 110
c            if not return to stmt. 15 and continue
           go to 15
c<<<<<<<<<<<<<<<<<<< END OF SORTING LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

100      continue
110   nfreq = ifreq - 1
      return
      end
