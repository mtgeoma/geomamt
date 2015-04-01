      subroutine ediformatfin( ind, permin, permax, nf )
c------------------------------------------------------------------------
c subroutine to return the number of frequencies between permin and permax
c in an EDI-format file
c------------------------------------------------------------------------

      implicit none
      
c------------------------------------------------------------------------

      include 'size.inc'
      
c------------------------------------------------------------------------

      character line*80
     
      real permin, permax, freqs(MAXDAT)
     
      integer nperdat, nf, i, ind
           
      logical DEBUG, INC
 
      common /cmndbg/ DEBUG
    
c------------------------------------------------------------------------
c...get number of periods in file (from entry after NFREQ=)
      read(ind,'(a)') line
      do while(index(line,'NFREQ=').eq.0)
        read(ind,'(a)') line
      enddo
      read(line(index(line,'=')+1:),*) nperdat    
      
      if( DEBUG ) write(*,*)'ediformatfin: No. periods in file: ', 
     &                      nperdat
       
       
      if( nperdat.gt.MAXDAT ) then
        write(*,*) 'Too many periods in file - greater '//
     &             'than MAXDAT'
        stop
      endif


c...look for >FREQ section
      do while(index(line,'>FREQ ').eq.0)
        read(ind,'(a)',end=10) line
      enddo
      if( index(line,'ORDER=INC').ne.0 ) then
        INC = .TRUE.
      else
        INC = .FALSE.
      endif

c...read in freqs
      if( DEBUG ) write(*,*)'ediformatfin: MTPARAMS edi file to be read'
      read(ind,*) (freqs(i),i=1,nperdat)
      goto 20

c...SPECTRA file. Have to get freqs off each SPECTRA section
10    if( DEBUG ) write(*,*)'ediformatfin: SPECTRA edi file to be read'
      rewind(unit=ind)
      do i = 1, nperdat

        read(ind,'(a)') line
        do while(index(line,'>SPECTRA ').eq.0)
          read(ind,'(a)') line
        enddo
        read(line(index(line,'FREQ=')+5:),*) freqs(i)
      enddo


c...convert to periods
20    do i = 1, nperdat
        freqs(i) = 1./freqs(i)
      enddo

c...check for number of periods between permin and permax
      nf = 0
      do i = 1, nperdat
        if(freqs(i).ge.permin .and. freqs(i).le.permax ) then
          nf = nf + 1
        endif
      enddo

      if( DEBUG ) write(*,*)'ediformatfin: nf =', nf  

100   return
      end
