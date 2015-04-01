      subroutine jformatfin( ind, permin, permax, nf )
c------------------------------------------------------------------------
c subroutine to return the number of frequencies between permin and permax
c in a J-format file
c------------------------------------------------------------------------

      implicit none
      
c------------------------------------------------------------------------

      include 'size.inc'
      
c------------------------------------------------------------------------

      character line*80
     
      real permin, permax, per
     
      integer nperdat, nf1, nf2, nf, idat, ind
           
      logical DEBUG

      common /cmndbg/ DEBUG

c------------------------------------------------------------------------
 
c...get number of periods in file (from entry after RXX or ZXX)
      read(ind,'(a)') line
      do while(index(line(2:3),'XX').eq.0)
        read(ind,'(a)') line
      enddo
      read(ind,*) nperdat    
      
      if( DEBUG ) write(*,*)'jformatfin: No. periods in file: ', nperdat
       
       
      if( nperdat.gt.MAXDAT ) then
        write(*,*) 'Too many periods in file - greater '//
     &             'than MAXDAT'
        stop
      endif
       
      nf1 = 0
      nf2 = 0
      do idat = 1, nperdat
        read(ind,*) per
        if( per.lt.0. ) per = -1./per

        if((per.ge.permin) .and. (per.le.permax)) nf = nf+1
      enddo

      if( DEBUG ) write(*,*)'jformatfin: nf =', nf  

      return
      end
