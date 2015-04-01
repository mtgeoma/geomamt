c
c******************************************************
c
      subroutine rfhead(iounit,nch,ndmax,nd,nfreqmx,nwin,dr
     1   ,chid,orient,decl,stcor,irecl)

c       writes header on fm and fg files; header is 20 records
c       ( = 160 bytes for fg (i.e. gds) ; = 240 bytes for fm (mt)
c          files).  

      include '../include/nchmx.inc'
      integer maxHeaderLength
      parameter (maxHeaderLength = 20*(nchmx*8+4))

      integer nwin(ndmax),idl(20)
      real dr(ndmax),stcor(2),orient(2,*),decl
      character*1 ctemp(maxHeaderLength)
      character*6 chid(*)
      character*12 cr

      do 200 i = 1,ndmax
      nwin(i) = 0
      dr(i) = 0.0
200   continue

      do 30 i = 1,20
      i1 = (i-1)*irecl + 1
      i2 = i1 + irecl-1
      read(iounit,rec = i) (ctemp(j),j=i1,i2)
c      write(6,'(16a1)') (ctemp(j),j=i1,i2)
30    continue

      i1 = 1
      call movec(ctemp,cr,8,i1)
      read(cr,101) nch
101   format(3x,i5)

      call movec(ctemp,cr,8,i1)
      read(cr,102) nd
102   format(2x,i6)

      call movec(ctemp,cr,10,i1)
      read(cr,103) nfreqmx
103   format(4x,i6)

      i1 = i1 + 4
      do 10 i = 1,nd
      call movec(ctemp,cr,2,i1)
      read(cr,'(i2)') idl(i)
      call movec(ctemp,cr,6,i1)
      read(cr,'(i6)') nwin(idl(i))
10    continue

      i1 = i1 + 4
      do i = 1,nd
         call movec(ctemp,cr,12,i1)
         read(cr,'(1pe12.6)') dr(idl(i))
      enddo

      i1 = i1 + 8
      do i = 1,nch
         call movec(ctemp,cr,12,i1)
         read(cr,'(a6,f6.1)') chid(i),orient(1,i)
         call movec(ctemp,cr,8,i1)
         read(cr,'(f8.2)') orient(2,i)
      enddo

      i1 = i1 + 4
      call movec(ctemp,cr,10,i1)
      read(cr,'(f10.5)') stcor(1)

      call movec(ctemp,cr,10,i1)
      read(cr,'(f10.5)') stcor(2)

      i1 = i1 + 4
      call movec(ctemp,cr,8,i1)
      read(cr,'(f8.4)') decl

      return
      end
