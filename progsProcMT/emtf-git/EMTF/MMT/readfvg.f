c_____________________________________________________________________
c
      subroutine readfvg(id,ifreq,ista,itape,nd,nsta,ntape,
     1  nch,xx,ix,nsets,ncht,maxget)
 
c       reads from FC files and returns data for frequency ifreq,
c       decimation level id, tape itape, station ista packed in array ix;
c       uses directory array irecd and input units array iounits to find
c       data; FIRST record after header contains scaling factor
c       This version unpacks the data into complex FCs
 
      include 'iosize.inc'

      integer ipack(nchmx),ix(nsta,*)
      complex xx(ncht,*),xt(nchmx)
      real xscale(nchmx)

      irec = irecd(id,ifreq,ista,itape)
      if(irec.gt.0) then
         read (inunit(ista,itape),rec = irec) id1,ifreq1,nsets
         nsets = min(nsets,maxget)
         if(lpack) then
            irec = irec+1
            read(inunit(ista,itape),rec = irec) (xscale(k),k=1,nch)

            do 15 i = 1,nsets
               irec=irec+1
               read (inunit(ista,itape),rec = irec)
     &                        ix(ista,i),(ipack(k),k=1,nch)
               call unpk(ipack,nch,xt)
               do 10 k = 1,nch
10                xx(k,i) = xt(k)*xscale(k)/1000.
15             continue
         else
            do 20 i = 1,nsets
               irec = irec + 1
               read(inunit(ista,itape),rec = irec ) 
     &                        ix(ista,i),(xx(k,i),k=1,nch)
20             continue
         end if
      else
         xx(1,1) = 0
         nsets = 0
      end if
      return
      end
c___________________________________________________________________
      subroutine unpk(ipack,n,x)
      integer ipack(n)
      complex x(n)

      do 15 j=1,n
         it1=mod(ipack(j),65536)
         if(it1.lt.0) it1=it1+65536
         it2=(ipack(j)-it1)/65536
         it1=it1-32768
         k1=mod(it1,8)
         if(k1.lt.0) k1=k1+8
         k2=mod(it2,8)
         if(k2.lt.0) k2=k2+8
         xt1=float(it1-k1)/32768.+.000122
         xt2=float(it2-k2)/32768.+.000122
         do 5 ii=1,k1
5           xt1=xt1*10.
         do 10 ii=1,k2
10          xt2=xt2*10.
         x(j)=cmplx(xt1,xt2)
15       continue
      return
      end
c______________________________________________________________________
c
      subroutine rfhead(iounit,nch,ndmax,nd,nfreqmx,nwin,dr
     1   ,chid,orient,decl,stcor,irecl)
 
c       reads header on fm and fg files; header is 20 records
c       ( = 160 bytes for fg (i.e. gds) ; = 240 bytes for fm (mt)
c          files).  
 
      include 'nstamx.inc'
      integer maxHeadLength
      parameter (maxHeadLength = 20*(8*nchmx+4))

      integer nwin(ndmax),idl(20)
      real dr(ndmax),stcor(2),orient(2,*),decl
      character*1 ctemp(maxHeadLength)
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
      do i = 1,nd
         call movec(ctemp,cr,2,i1)
         read(cr,'(i2)') idl(i)
         call movec(ctemp,cr,6,i1)
         read(cr,'(i6)') nwin(idl(i))
      enddo

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
c______________________________________________________________________
c
      subroutine movec(c1,c2,nmv,i1)
      character*1 c1(*)
      character*120 c2

      do 5 i = 1,nmv
      j = i1 - 1 + i
      c2(i:i)=c1(j)
5     continue

      i1 = i1 + nmv

      return
      end
