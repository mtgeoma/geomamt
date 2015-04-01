c______________________________________________________________________
c
      subroutine mkfdir(nsta,ntape)
 
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
         irec = nhdrec + 1
         do 50 i = 1,ndmax
         do 50 j = 1,nfreqmax
            if (lfop)then
               open(unit=iounit,file=cfilein(ista,itape),
     &         form = 'unformatted',access='direct',recl=iorecl(ista))
            end if
            read(iounit,rec=irec,iostat=ioerr) id,ifreq,nsets
            if(ioerr.ne.0) go to 90
            if(lfop) close(iounit)
            if((id.gt.ndmax).or.(ifreq.gt.nfreqmax).or.(id.lt.0).or.
     1         (ifreq.lt.0)) then
               print*,'error; id = ',id,'ndmax = ',ndmax,'ifreq = ',
     1         ifreq,'nfreqmax = ',nfreqmax
               ifn = (i-1)*nfreqmax+j
               print*,'ista,frequency #',ista,ifn
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
