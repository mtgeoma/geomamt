C ME    REMOVE SCRATCH FILE
c      subroutine outinit(iouout,ioscr,nch,cfout,lpack)
C ME    REMOVE SCRATCH FILE
      subroutine outinit(nch,lpack)
cnew      character*80 cfout

      include 'iounits.inc'
      include '../include/four_byte.inc'
      logical lpack
c       initialize for outputing packed record


C ME    REMOVE SCRATCH FILE
c       open scratch file for temporary output of FCs
c      irlscr = 4*(2*nch+1)
c      if(l_4byte) irlscr = irlscr/4
c      open(file='f_scratch',unit=ioscr,form='unformatted',
c     &    access='direct',recl=irlscr)
C ME    REMOVE SCRATCH FILE

c        open file for output of ordered FCs
      if(lpack) then
         irecl = 4*(nch+1)
      else
         irecl = 4*(2*nch+1)
      end if
      if(l_4byte) irecl = irecl/4

      open(unit=out_unit,file=cfout,access='direct',form='unformatted',
     &  recl=irecl)
      return
      end
c________________________________________________________________
c

C     ME        REMOVE SCRATCH FILE
c      subroutine mkseq(iouin,iouout,ifdir,iset,nsets,nch,
c     1  nd,cwrk,cout,nfmax,idoff,pspecl1,nchmx,nwmx,lpack,iwrk,ixs)
       subroutine mkseq(ifdir,iset,nsets,nch,
     1  nd,cwrk,cout,nfmax,idoff,pspecl1,nchmx,nwmx,lpack,iwrk,ixs,
     2     nsmax,x_scr, iuse_scr) 
        real x_scr(nsmax*nfmax,nch,2)
        integer iuse_scr(nsmax*nfmax)
C     ME        REMOVE SCRATCH FILE


c    overly complex, because it can treat the case of missing frequencies
c       makes sequential file ordered by frequency (i.e. all fourier
c       coefficients for a given frequency are adjacent, ordered by
c       set number)
      complex cwrk(nch,nsets),cout(nch,nsets)
      integer iset(4,nsets),ifdir(3,nsets),iwrk(nch,nsets),ixs(nsets)
      real pspecl1(nchmx,nwmx,*)
      logical lpack  !lhead,

      ireco = 20 
      do 50 id=idoff + 1,idoff + nd
c       set up; lrec  - last record for a set; ndr  - no. of sets for
c       decimation level id
         lrec=0
         ndr=0
         do 10 is=1,nsets
            irec=lrec+1
            lrec=lrec+ifdir(2,is)
            if(ifdir(1,is).ne.id) go to 10
            ndr=ndr+1

c       iset - for each set for decimation level id, iset contains:
c               1. set no.
c               2. frequency no.
c               3. last record read for this set - stored in cwrk
c               4. final record to read
            iset(4,ndr)=lrec
            iset(3,ndr)=irec
            iset(1,ndr)=ifdir(3,is)
C     ME        REMOVE SCRATCH FILE
c            read(iouin,rec=irec) iset(2,ndr),(cwrk(k,ndr),k=1,nch)
        iset(2,ndr) = iuse_scr(irec)     ! SCRATCH FILE 02.01.98
        do k = 1,nch                                 ! SCRATCH FILE 02.01.98
        cwrk(k,ndr) = cmplx(x_scr(irec,k,1),x_scr(irec,k,2))     ! SCRATCH FILE 02.01.98
        enddo                                                        ! SCRATCH FILE 02.01.98

10          continue

c       sort file and output
        do 40 j=1,nfmax
           nf=0
           do 30 i=1,ndr
              if(iset(2,i).eq.j) then
                if(iset(3,i).le.iset(4,i)) then

                   nf=nf+1
                   do 20 k = 1,nch
20                 cout(k,nf) = cwrk(k,i)*1000./pspecl1(k,j,id-idoff)
                   ixs(nf) = iset(1,i)

                   if(iset(3,i).lt.iset(4,i)) then
c                    not done; read another record from set 
                       iset(3,i)=iset(3,i)+1
C     ME        REMOVE SCRATCH  FILE
c                       read(iouin,rec=iset(3,i)) iset(2,i),
c     &                               (cwrk(k,i),k=1,nch)
        i_scr_rec = iset(3,i)
        iset(2,i) = iuse_scr(i_scr_rec) 
        do k = 1,nch                  
           cwrk(k,i) = cmplx(x_scr(i_scr_rec,k,1),x_scr(i_scr_rec,k,2))
        enddo
C     ME        REMOVE SCRATCH  FILE

                    end if

                 end if
              end if
30            continue
           call pack(nch,nf,cout,ixs,j,id,ireco,
     &                      pspecl1(1,j,id-idoff),lpack,iwrk)
40         continue
50       continue
      return
      end
c
c______________________________________________________________________
c
      subroutine pack(nch,nf,cout,ixs,ifreq,id,irec,pspecl1,
     &   lpack,ix)

      include 'iounits.inc'
c    outputs both real and imaginary parts of complex fourier coefficients
c           packed into one four byte integer 

      parameter (minint = -2147483647)
      real  cout(2,nch,nf),pspecl1(nch)
      integer ix(nch,nf),nch,nf,id,irec,ixs(nf)
      logical lpack

      irec = irec + 1
      write(out_unit,rec=irec) id,ifreq,nf

      if(lpack) then
         irec = irec + 1
         write(out_unit,rec=irec) pspecl1
         do 20 i=1,nf
         do 15 j=1,nch
            xt1=cout(1,j,i)
            xt2=cout(2,j,i)
               do 5 k=1,8
               if(abs(xt1).lt.(1.0)) go to 6
5              xt1=xt1/10.0
            ix(j,i)=minint
             print*,'error in packing'
            go to 15
6           k1=k-1
            l1=nint(4096*xt1-.5)
               do 10 k=1,8
               if(abs(xt2).lt.(1.0)) go to 11
10             xt2=xt2/10.0
            ix(j,i)=minint
             print*,'error in packing'
            go to 15
11          k2=k-1
            l2=nint(4096*xt2-.5)
            ix(j,i)=65536*(8*l2+k2)+32768+(8*l1+k1)
15          continue
20       continue

         do 25 i = 1,nf 
            irec = irec+1
            write(out_unit,rec=irec) ixs(i),(ix(j,i),j=1,nch)
25          continue
      else
         do 30 i = 1,nf
            irec = irec+1
            write(out_unit,rec=irec) 
     &           ixs(i),((cout(k,j,i),k=1,2),j=1,nch)
30          continue
      end if
      return
      end
c______________________________________________________________________
c
      subroutine wfhead(nch,nd,nfreqmx,nwin,dr,idl
     1   ,chid,orient,decl,stcor,irlo)

      include 'iounits.inc'
      include '../include/nchmx.inc'
      integer maxHeaderLength 
      parameter (maxHeaderLength = 20*(nchmx*8+4))

      integer nwin(nd),idl(nd)
      real dr(nd),stcor(2),orient(2,nch),decl
      character*1 ctemp(maxHeaderLength)
      character*6 chid(nch)
      character*12 cr

ccc   ctemp of this size should allow for 100 channels + 20 decimation levels
      data ctemp/maxHeaderLength*' '/

      if(nd>3*nch) then
         write(6,*) "nd>3*nch cannot be properly writing by wfhead"
         stop
      endif
      if(nd>99) then
         write(6,*) "nd>99 cannot be properly writing by wfhead"
         stop
      endif

      i1 = 1
      write(cr,101) nch
101   format('nch',i5)

      call cmove(ctemp,cr,8,i1)

      write(cr,102) nd
102   format('nd',i6)

      call cmove(ctemp,cr,8,i1)

      write(cr,103) nfreqmx 
103   format('nfmx',i6)

      call cmove(ctemp,cr,10,i1)

      call cmove(ctemp,'nwin',4,i1)
      do i = 1,nd
         write(cr,'(i2i6)') idl(i),nwin(i)
         call cmove(ctemp,cr,8,i1)
      enddo

      call cmove(ctemp,'dr  ',4,i1)
      do i = 1,nd
         write(cr,'(1pe12.6)') dr(i)
         call cmove(ctemp,cr,12,i1)
      enddo

      call cmove(ctemp,'ch,orien',8,i1)
      do i = 1,nch
         write(cr,'(a6,f6.1)') chid(i),orient(1,i)
         call cmove(ctemp,cr,12,i1)
         write(cr,'(f8.2)') orient(2,i)
         call cmove(ctemp,cr,8,i1)
      enddo

      call cmove(ctemp,'scor',4,i1)
      write(cr,'(f10.5)') stcor(1)

      call cmove(ctemp,cr,10,i1)
      write(cr,'(f10.5)') stcor(2)

      call cmove(ctemp,cr,10,i1)

      call cmove(ctemp,'decl',4,i1)
      write(cr,'(f8.4)') decl

      call cmove(ctemp,cr,8,i1)

      do i = 1,20
         i1 = (i-1)*irlo + 1
         i2 = i1 + irlo-1
         write(out_unit,rec=i) (ctemp(j),j=i1,i2)
      enddo

      return
      end
