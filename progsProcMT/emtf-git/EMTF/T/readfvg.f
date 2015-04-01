c_____________________________________________________________________
c
      subroutine readfvg(id,ifreq,ista,itape,nd,nsta,ntape,
     1  nch,xx,ix,nsets,ncht)
 
c       reads from fg and fm files and returns data for frequency ifreq,
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
