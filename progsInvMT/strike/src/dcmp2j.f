c********************************************************************
c
c   dcmp2j
c
c   This program reads strike output file (.dcmp) 
c   and converts it to J-format  (******g.dat)
c    
c********************************************************************
c
c  I/O
c
c  Unit 10: Input file  (xxxxxx.dcmp)
c  Unit 12: Output file (xxxxxxg.dat)
c
c------------------------------------------------------------------------

      implicit none

c------------------------------------------------------------------------
 
      include 'size.inc'
      include 'version.inc'
      include 'ctrlblk.inc'

c------------------------------------------------------------------------

      real period(MAXF), per
      real rhoa(MAXF),phaa(MAXF),rhob(MAXF),phab(MAXF)
      real rhoamax(MAXF),rhoamin(MAXF),phaamax(MAXF),phaamin(MAXF)
      real rhobmax(MAXF),rhobmin(MAXF),phabmax(MAXF),phabmin(MAXF)
      real xave(5), x(5,MAXF)
      real azimuth, anis, azimdat

      integer i, j, k, l, ifact
      integer npoints, n_an1, n_an2, ndatum
      integer ltitles(11)
      integer lnblnk

      character filename*50, outfile*50, station*40
      character*24 date_string
      character titles(11)*20
      character*80 line

      character*3 type(4)/'RXY','RYX','RRG','RRD'/
      character*1 chk/'Y'/, chkt/'Y'/, quadflip/'N'/

      LOGICAL QANIS, WAZIM

c-------------------------------------------------------------------------------
c   initialize

      call fdate( date_string )

      ifact=2
      ndatum=1

      titles(1)  = 'regional azimuth'
      titles(2)  = 'shear angle'
      titles(3)  = 'channelling angle'
      titles(4)  = 'twist angle'
      titles(5)  = 'app rho a'
      titles(6)  = 'app rho b'
      titles(7)  = 'imped phase a'
      titles(8)  = 'imped phase b'
      titles(9)  = 'error'
      titles(10) = 'skew'
      titles(11) = 'anisotropy'

      do i = 1, 11
        ltitles(i)  = lnblnk(titles(i))
      enddo

c-------------------------------------------------------------------------------


      write(*,*)' '
      write(*,*)'                Welcome to dcmp2j'
      write(*,*)'                  version: ', ver
      write(*,*)'        Todays date: ', date_string
      write(*,*)' '
      write(*,*)'               Version limits:'
      write(*,*)'      Max. no. freqs:       ', MAXF
      write(*,*)' '
      write(*,*)'     This program reads output (.dcmp) from STRIKE: '
      write(*,*)'     Places in standard J-format >>   "______g.dat'
      write(*,*)' '

c-------------------------------------------------------------------------------

      iprint = 0
      call iin('Give print level (0/1/2/5/10)', iprint )

c-------------------------------------------------------------------------------

c...get station name
      station='tbt101'
      call tin('Station name?', station)

c...open input and output files
      filename= station(1:lnblnk(station))//'.dcmp'
      outfile = station(1:lnblnk(station))//'g.dat'
      open(unit=10,file=filename,status='old' )
      open(unit=12,file=outfile)

      if( iprint.ge.0 ) then
        write(*,*)'Opened input  file >',filename(1:lnblnk(filename))
        write(*,*)'Opened output file >',outfile(1:lnblnk(outfile))
      endif

c...get azimuth of original data
      read(10,'(a)') line
      do while( line(1:8).ne.'>AZIMUTH' )
        read(10,'(a)',end=1) line
      enddo
      WAZIM = .TRUE.
      read(line(13:),*) azimdat
      if(iprint.ge.1) write(*,*)'Azimuth of original data =', azimdat
      goto 2
1     WAZIM = .FALSE.
      write(*,*)'Keyword AZIMUTH not found in .dcmp file -'//
     &          ' assumed to be zero'
2     continue

c...request whether to perform anisotropy correction
      QANIS = .TRUE.
      call lin('Correct anisotropy', QANIS )

c...get average values of azimuth, shear, chann & twist
      read(10,'(a)') line
      do k=1,4
        do while(index(line,titles(k)(1:ltitles(k))).eq.0)
          read(10,'(a)',end=9001) line
        enddo
        if( iprint.ge.2 ) then
          write(*,*) 'Block ',titles(k),' found'
          write(*,*) 'Line >', line(1:lnblnk(line))
        endif
        read(line,*) npoints
        if( iprint.ge.2 ) then
          write(*,*) 'npoints =', npoints
        endif
        xave(k)=0.
        do l=1,npoints
          read(10,*) period(l), x(k,l)
          xave(k)=xave(k) + x(k,l)
        enddo
        xave(k)=xave(k)/npoints
      enddo

c...get average value of error
      k = 9
      do while(index(line,titles(9)(1:ltitles(9))).eq.0)
        read(10,'(a)',end=9001) line
      enddo
      read(line,*) npoints
      xave(5)=0.
      do l=1,npoints
        read(10,*) period(l), x(5,l)
        xave(5)=xave(5) + x(5,l)
      enddo
      xave(5)=xave(5)/npoints
      rewind(10)

c...write out comment block
      write(12,'(a)') '# Written by dcmp2j: input file >'//
     &                   filename(1:lnblnk(filename))
      write(12,'(a)') '# date: '//date_string
      write(12,'(a)') '#'
      write(12,'(a)') '#'
      write(12,'(a,f5.1)') '# azimuth of original data:', azimdat
      write(12,'(a)') '#'
      write(12,'(a)') '# The azimuth listed under AZIMUTH is the'//
     &                  ' average value in above azimuth coordinates'
      write(12,'(a)') '#'
      write(12,'(a)') '# average distortion parameters:'
      write(12,'(a,f5.1)') '# azimuth:', xave(1)
      write(12,'(a,f5.1)') '# shear  :', xave(2)
      write(12,'(a,f5.1)') '# twist  :', xave(4)
      write(12,'(a,f5.1)') '# chann  :', xave(3)
      write(12,'(a,f5.1)') '# error  :', xave(5)
      write(12,'(a)') '#'
      write(12,'(a)') '#     Period       Azim    Shear'//
     &                '   Chann   Twist     Error'
      do l = 1, npoints
        write(12,'(a,f12.4,5f8.1)') '#',period(l),(x(k,l),k=1,5)
      enddo
      write(12,'(a)') '#'


c...copy comment and info blocks to output file
      read(10,'(a)') line
      do while(line(1:1).eq.'#')
        write(12,'(a)') line(1:lnblnk(line))
        read(10,'(a)') line
      enddo
      backspace(unit=10)


c...copy info blocks, changing azimuth info to that of average decomp direction
      read(10,'(a)') line
      do while(line(1:1).eq.'>')
        if(index(line,'>AZIMUTH').ne.0 )then
          write(line(13:),*) xave(1)
          write(*,*)'Written AZIMUTH of', xave(1)
        endif
        write(12,'(a)') line(1:lnblnk(line))
        read(10,'(a)') line
      enddo

c...need to write AZIMUTH if not written to this time
      if( .not.WAZIM ) then
        line = '>AZIMUTH'
        write(line(13:),*) xave(1)
        write(12,'(a)') line(1:lnblnk(line))
      endif

      write(12,555) outfile,(xave(k),k=1,4)
555   format(A11,3x,'AZIM',F8.2,' SHEAR',F8.2,'  CHANN',
     >         F8.2,'  TWIST',F8.2)
      write(*,*)' '
      write(*,*)' Azimuth of X-axis is :>>  ',xave(1)

c...app rho a
      rewind( 10 )
      do while( index(line,'app rho a').eq.0 )
        read(10,'(a)') line
      enddo
      read(line,*)npoints
      do j=1,npoints
        read(10,*)period(j),rhoa(j),rhoamin(j), rhoamax(j)
        if( rhoamin(j).eq.0.) rhoamin(j) = -999.
        if( rhoamax(j).eq.0.) rhoamax(j) = -999.
      enddo

c...imped phase a
      rewind( 10 )
      do while( index(line,'imped phase a').eq.0 )
        read(10,'(a)') line
      enddo
      read(line,*)npoints
      do j=1,npoints
        read(10,*)period(j),phaa(j),phaamin(j),phaamax(j)
        if( phaamin(j).eq.0.) phaamin(j) = -999.
        if( phaamax(j).eq.0.) phaamax(j) = -999.
      enddo

c...app rho b
      rewind( 10 )
      do while( index(line,'app rho b').eq.0 )
        read(10,'(a)') line
      enddo
      read(line,*)npoints
      do j=1,npoints
        read(10,*)period(j),rhob(j),rhobmin(j), rhobmax(j)
        if( rhobmin(j).eq.0.) rhobmin(j) = -999.
        if( rhobmax(j).eq.0.) rhobmax(j) = -999.
      enddo

c...imped phase b: move into III quadrant
      rewind( 10 )
      do while( index(line,'imped phase b').eq.0 )
        read(10,'(a)') line
      enddo
      read(line,*)npoints
      do j=1,npoints
        read(10,*)period(j),phab(j),phabmin(j),phabmax(j)
        phab(j) = phab(j) - 180.
        if( phabmin(j).eq.0.) then
          phabmin(j) = -999.
        else
          phabmin(j) = phabmin(j) - 180.
        endif
        if( phabmax(j).eq.0.) then
          phabmax(j) = -999.
        else
          phabmax(j) = phabmax(j) - 180.
        endif
      enddo


c...correct anisotropy
      if( QANIS ) then
        write(*,*) '  #           Period        Anisotropy'
        do j = 1, npoints
          per = period(j)
          if( per.lt.1. ) per = -1./per
          write(*,'(i5,f15.1,f15.3)') j, per, sqrt(rhoa(j)/rhob(j))
        enddo
100     call nin('Give period range to average (e.g. 1,5)')
        read(*,*) n_an1, n_an2
        if( n_an1.le.0 .or. n_an1.gt.n_an2 .or.
     &      n_an2.gt.npoints ) then
          write(*,*)'Inappropriate values entered: re-enter'
          goto 100
        endif
        anis = 0.
        do j = n_an1, n_an2
          anis = anis + sqrt(rhoa(j)/rhob(j))
        enddo
        anis = anis/(float(n_an2 - n_an1 + 1))
        write(*,*) 'Average anisotropy to be applied:', anis
        do j = 1, npoints
          rhoa(j) = rhoa(j)/anis
          if( rhoamin(j).ne.-999.)rhoamin(j) = rhoamin(j)/anis
          if( rhoamax(j).ne.-999.)rhoamax(j) = rhoamax(j)/anis
          rhob(j) = rhob(j)*anis
          if( rhobmin(j).ne.-999.)rhobmin(j) = rhobmin(j)*anis
          if( rhobmax(j).ne.-999.)rhobmax(j) = rhobmax(j)*anis
        enddo
      endif

c...write data out
      write(12,'(a)') 'RXY'
      write(12,*) npoints
      do i = 1, npoints
        write(12,2002) period(i), rhoa(i), phaa(i),
     &                 rhoamin(i), rhoamax(i),
     &                 phaamin(i), phaamax(i),
     &                 1., 1.
      enddo

      write(12,'(a)') 'RYX'
      write(12,*) npoints
      do i = 1, npoints
        write(12,2002) period(i), rhob(i), phab(i), 
     &                 rhobmin(i), rhobmax(i),
     &                 phabmin(i), phabmax(i),
     &                 1., 1.
      enddo

      stop

9001  write(*,*)'EOF found looking for block keyword <', titles(k),'>'
      stop

9002  write(*,*)'EOF found looking for AZIMUTH'
      stop



2002  format( g11.5,2x,g11.5,1x,f6.1,1x,2(g10.4,1x),1x,
     &                    2(f6.1,1x),2(f4.2,1x) )

      end
   
