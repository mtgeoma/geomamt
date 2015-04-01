      subroutine fop(nsta,nch,ntape)
      include 'iosize.inc'
      integer ntape(*),nch(*)
c   open fourier coefficient files as direct access unformatted files       
      do 2 i=1,nsta
         do 2 j=1,ntape(i)
            open(unit=inunit(i,j),file=cfilein(i,j),form='unformatted'
     1                ,access='direct',recl=iorecl(i),status = 'old')
c            write(*,*) 'unit   ',inunit(i,j),'   file   ',
c     &           cfilein(i,j),'open'
2           continue
      return
      end 
