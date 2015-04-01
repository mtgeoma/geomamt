c______________________________________________________________________
c
      subroutine mk_list(nd,id,nsta,ntape,nch,nsets,isets,cwrk)
 
ccc   reads in all files in cf_array file and makes a list of all
ccc   available set numbers for decimation level id
ccc   outputs # of sets for each station in nsets, and list of actual
ccc   set numbers in array isets

c     If lfop is .true. open and close files for each call to readfvg
c      use this way when number of files is large

      include 'iosize.inc'

      integer ntape(*),isets(nsta,*),nch(nsta),nsets(nsta)
      complex cwrk(*)

      data maxint/2147483647/

      ib = 1
c********* read loop for nsta stations
      do ista = 1,nsta
         i0 = 1
         do itape = 1,ntape(ista)
            if(lfop) then
c              need to open file before reading
               open(unit=inunit(ista,itape),
     &             file=cfilein(ista,itape),form='unformatted',
     &             recl=iorecl(ista),access='direct')
             end if
ccc          maxget is set to insure that the total number of
ccc          points doesn't exceed nsetmx (set in iosize.h)
             maxget = ntfmax - i0 + 1
             write(*,*) 'id,ib,ista,itape,nd,nsta'
             write(*,*) id,ib,ista,itape,nd,nsta
             call readfvg(id,ib,ista,itape,nd,nsta,ntape,
     &         nch(ista),cwrk,isets(1,i0),nf,nch(ista),maxget)
c                 write(*,*) 'NF,ISTA,I0,MAXGET',nf,ista,i0,maxget
             if(lfop) close(inunit(ista,itape))
             i0 = i0 + nf
          enddo ! itape = 1,ntape(ista)
          nsets(ista) = i0 - 1
       enddo   ! ista = 1,nsta
       return
       end
ccc_____________________________________________________________________
ccc
      subroutine wrt_sets(nsta,id,sta,nsets,isets,cfile)
      character*80 cfile
      character*3 sta(nsta)
      include 'nstamx.inc'
      integer nsta,nsets(nsta),isets(nsta,*),fid,id
      integer iseg(2,nstamx)
      fid=55
ccc      open(unit= fid,form='unformatted',file=cfile)
      open(unit= fid,form='formatted',file=cfile)
      write(fid,*) nsta,id
      write(fid,'(10(a3,1x))') sta
      write(fid,*) nsets
      do ista = 1,nsta
         nseg = 1
         l = 1
         iseg(1,1) = isets(ista,1)
         do k= 2,nsets(ista)
            if(isets(ista,k) .eq.iseg(1,nseg)+l) then
ccc            data sets are continuous, add one to l
               l = l + 1
            else
ccc            gap in data sets; terminate old segment; start new segment
               iseg(2,nseg) = iseg(1,nseg)+l-1
               nseg = nseg + 1
               iseg(1,nseg) = isets(ista,k)
               l = 1
            endif
          enddo
          iseg(2,nseg) = isets(ista,nsets(ista))
          call wrt_1(fid,nseg,iseg)
       enddo
       close(fid) 
       return
       end
ccc_____________________________________________________________________
ccc
       subroutine wrt_1(fid,nseg,iseg)
       integer fid,nseg,iseg(2,nseg)

       write(fid,'(i6)') nseg
       do k=1,nseg
          write(fid,'(2i6)') (iseg(l,k),l=1,2)
       enddo
       return
       end
