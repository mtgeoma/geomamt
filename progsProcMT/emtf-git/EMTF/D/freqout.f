c_____________________________________________________________
c
C     ME	remove SCRATCH FILE
c      subroutine freqout(iounit,nch,nf,x,irec)   
      subroutine freqout(nch,nf,x,irec,nsmax,nfusemx,
     &  x_scr,iuse_scr)   
C     ME	remove SCRATCH FILE

c    new version doesn't do any packing; just outputs
c      FCs to scratch file
 
      real  x(nch,2,nf)
 
C   SCRATCH FILE 02.01.98
        real x_scr(nsmax*nfusemx,nch,2)
        integer iuse_scr(nsmax*nfusemx)
C   SCRATCH FILE 02.01.98
      do 30 i = 1,nf
      irec = irec+1
C   SCRATCH FILE 02.01.98
c      write(iounit,rec=irec) i,(x(j,1,i),x(j,2,i),j=1,nch)
        iuse_scr(irec) = i  !   SCRATCH FILE 02.01.98
        do j = 1,nch                             !   SCRATCH FILE 02.01.98
                x_scr(irec,j,1) = x(j,1,i)  !   SCRATCH FILE 02.01.98
                x_scr(irec,j,2) = x(j,2,i)  !   SCRATCH FILE 02.01.98
        enddo                                    !   SCRATCH FILE 02.01.98
C   SCRATCH FILE 02.01.98
30    continue
      return
      end
c_____________________________________________________________
c
      subroutine l1spec(nch,x,nf,pspecl1,nchmx)

c    new version doesn't do any packing; just outputs
c      FCs to scratch file
      real pspecl1(nchmx,*),x(nch,2,nf),temp
 
      do 30 i = 1,nf
      do 30 j = 1,nch
         temp = (x(j,1,i)**2.+x(j,2,i)**2.)
         if(temp.gt.0) then
            temp = sqrt(temp)
         else
c            print*,'temp = ',temp
c            print*,'sqrt(temp) = ',sqrt(temp)
            temp = 0.
         end if
         pspecl1(j,i) = pspecl1(j,i)+temp
30       continue
      return
      end
