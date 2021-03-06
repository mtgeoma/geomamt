c______________________________________________________________________
c
      subroutine mkfdir(nsta,ntape)
 
ccc     ARRAY VERSION
c       makes directory of FC files so that data from any given
c       frequency can be accessed directly; assumes direct access conection of
c       units in array iounits; header of length nhdrec records is assumed;
c       if lpack = .true. then an extra record (for variable scaling)
c        is assumed for each frequency block

      include 'iosize.inc'
      integer ntape(nsta)
          
      do 100 ista=1,nsta
         do 100 itape=1,ntape(ista)
            do 5 i=1,ndmax
            do 5 j=1,nfreqmax
5           irecd(i,j,ista,itape) = 0
         iounit = inunit(ista,itape)
c         write(*,*) 'ista,itape,iounit',ista,itape,iounit
         irec = nhdrec + 1
         do 50 i = 1,ndmax
         do 50 j = 1,nfreqmax
            if (lfop)then
               open(unit=iounit,file=cfilein(ista,itape),
     &         form = 'unformatted',access='direct',recl=iorecl(ista))
            end if
            read(iounit,rec=irec,iostat = ioerr) id,ifreq,nsets
            if(ioerr .ne. 0 ) go to 90
            if(lfop) close(iounit)
            if((id.gt.ndmax).or.(ifreq.gt.nfreqmax).or.(id.lt.0).or.
     1         (ifreq.lt.0)) then
               print*,'error; id = ',id,'ndmax = ',ndmax,'ifreq = ',
     1         ifreq,'nfreqmax = ',nfreqmax
               ifn = (i-1)*nfreqmax+j
               print*,'ista,frequency #,nsets',ista,ifn,nsets
  
            else
               if(nsets.gt.0) irecd(id,ifreq,ista,itape) = irec
            end if
            if(lpack) then
               irec=irec+nsets+2
            else
               irec=irec+nsets+1
            end if
50          continue  

c     after hitting eof, close file; note: file must be reopened before reading
90       close(iounit)
100   continue
      return
      end
c______________________________________________________________________
c
      subroutine mkrec(nd,id,iband,nsta,ntape,nch,ncht,isuse,
     &  nfreq,xx,llx,lx,ixs,ixf,cwrk,iwrk)
 
ccc     ARRAY VERSION
c  subroutine makes complex FC array records (stored in array xx) for ntape
c   tapes of data, nsta stations; record includes frequencies iband(1)
c   through iband(2) at decimation level id; requires the directory array
c   irecd (in commonblock ioblock) be constructed by a call to
c      routine mkfdir bdfore first call to mkrec;
c    other inputs: ncht = total channels of data
c     isuse: Only sets at decimation level id with numbers between
c         isuse(1,id) and isuse(2,id) will be output


ccc    IN THIS VERSION: llx is a logical paramter which controls
ccc      how the routine handles sets for which some stations do
cccc     not have data      
c    IF llx = .true.  the routine returns a logical array lx,
c       the elements of which are true or false depending on whether
c       good data is present for a given station at a given frequency;
c    In this case storage for the logical array lx of total dimension
c     nsta*nfreq must be provided for in the calling program;
c    if llx = .false. lx only has to be of dimension nsta
c      in this case only sets for which all data are present and good
c      will be returned; standard usage for routine TF processing
ccc    NOTE SLIGHT DIFFERENCE FROM mkrec routine used for standard TF
ccc    processing

c     if lfull = .true. return set number (in ixs) and frequency band
c       number (in ixf); otherwise just FCs

c     If lfop is .true. open and close files for each call to readfvg
c      use this way when number of files is large

      include 'iosize.inc'

      logical lx(nsta,*),lset(100),luse,llx
      integer ntape(*),
     1  iset(100),lstset(100),iband(2 ),ncht0(100),iwrk(nsta,*),
     2  isuse(2,nd),nch(nsta)
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
c      write(*,*) 'IBAND',iband(1),iband(2)
c      write(*,*) 'ND,ID',nd,id
c      write(*,*) 'ISUSE',(isuse(k,id),k=1,2)
c      write(*,*) 'NCH',nch
c      write(*,*) 'NTAPE',(ntape(k),k=1,nsta)
c      write(*,*) 'LFOP',lfop
c      write(*,*) 'inunit = ',(inunit(l,1),l=1,6)


      do 100 ifb=iband(1),iband(2)
         nmax=0
c********* read loop for nsta stations
          do ista = 1,nsta
             i0=1
             do itape = 1,ntape(ista)
ccc        maxget is set to insure that the total number of
ccc        points doesn't exceed nsetmx (set in iosize.h)
                maxget = nsetmx - i0 + 1
                if (lfop) then
c                 need to open file before reading
                   open(unit=inunit(ista,itape),
     &             file=cfilein(ista,itape),form='unformatted',
     &             recl=iorecl(ista),access='direct')
                end if
c                write(*,*) 'ISTA,IORECL',ista,iorecl(ista)
                call readfvg(id,ifb,ista,itape,nd,nsta,ntape,
     &            nch(ista),cwrk(ncht0(ista)+1,i0),
     &            iwrk(1,i0),nf,ncht,maxget)
c                 write(*,*) 'NF,ISTA,I0,MAXGET',nf,ista,i0,maxget
                if(lfop) close(inunit(ista,itape))
                i0 = i0 + nf
             enddo ! itape = 1,ntape(ista)
             nf = i0 - 1
c               nmax is maximum number of sets over all stations
             if (nf .gt. nmax) nmax = nf
c                 mark the end of set number array 
             iwrk(ista,nf+1)=maxint
c                initialize iset, lstset
             do k = 1,nf
                if(abs(iwrk(ista,k)).ge.isuse(1,id)) go to 8
             enddo
             k = nf+ 1
8            continue

             iset(ista)=abs(iwrk(ista,k))
             lset(ista) = (iwrk(ista,k).ge.0).or.ljunk(ista)
             lstset(ista)=k
          enddo   ! ista = 1,nsta
c          write(*,*) 'Starting Set numbers',iset(1),iset(2)
c*********  data is now read in for frequency ifb, and stored in array cwrk
c           Now: organize into records with all set numbers the same for
c               all stations
 
c>>>>>>>>>>>>>>>>>>>>>> START OF SORTING LOOP >>>>>>>>>>>>>>>>>>>>>>>>>>>
15         continue
c            find smallest set number
           nxtset=iset(1)
           do 20 ista=2,nsta
20            if(nxtset.gt.iset(ista)) nxtset=iset(ista)
 
c******** if nxtset .eq. maxint  then done with frequency ifb
           if((nxtset.eq.maxint).or.(nxtset.gt.isuse(2,id))) go to 100
 
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
110   continue
      nfreq = ifreq - 1
      return
      end
c___________________________________________________________________
      subroutine fop(nsta,nch,ntape)
      include 'iosize.inc'
      integer ntape(*),nch(*)
c   open fourier coefficient files as direct access unformatted files       
      do i=1,nsta
         do j=1,ntape(i)
            open(unit=inunit(i,j),file=cfilein(i,j),form='unformatted'
     1                             ,access='direct',recl=iorecl(i))
            read(inunit(i,j),rec=nhdrec+3) istart
ccc            write(36,*) i,j,cfilein(i,j),istart
         enddo
      enddo
      return
      end 
