c========================================================================
c
      program strike
c
c     Gary McNeice 
c     Geosystem Canada Ltd.
c     927 Raftsman Lane
c     Ottawa, Ontario, K1C 2V3
c     CANADA
c     gmcneice@netcom.ca
c
c     Alan Jones 
c     Dublin Institute for Advanced Studies
c     5 Merrion Square
c     Dublin 2
c     Ireland
c     alan@cp.dias.ie
c
c------------------------------------------------------------------------
c Conditions for the use of strike:
c
c 1) I accept the McNeice-Jones multi-site, multi-frequency magnetotelluric
c    tensor decomposition code, hereinafter called "strike", on a caveat
c    emptor basis
c
c 2) I will use strike for non-profit purposes only
c
c 3) I will not accept any payment for use of strike
c
c 4) I will not give strike to any other person
c
c 5) I will inform the authors of strike, Gary McNeice and Alan Jones, of any
c    coding errors that I find
c
c 6) I will inform the authors of strike, Gary McNeice and Alan Jones, of any
c    improvements and additions that I make
c
c 7) I will acknowledge the use of strike in any publications and presentations
c    I write or give in which I present results based on strike-based
c    decomposition or decomposed data
c
c 8) I understand that the current version of strike, v4.0, uses NAG routines
c    that I have to provide myself
c
c------------------------------------------------------------------------
c
c
c  v3.0: -reads in J-format (.dat) files and computes realizations internally
c        -does NOT read in .g file
c        -scales impedances by 1/sqrt(mu*omega)
c
c  v3.1: -choice of normalizing term in objective function (old=GAVSD2)
c
c  v3.2: -corrects format significance in imped*.dat output file
c        -writes out normalized rms error
c        -sets error floor to percentage of maximum value (EXPERT controllable)
c
c  v3.3: -writes out phase differences
c        -corrected minor bug when last file has no freqs
c        -increased filename size
c
c  v3.4: -minor mods for compilation with GNU g77 compiler
c
c  v4.0: -reads in .edi files OR .dat files
c
c  v4.1: -corrected azimuth error for EDI files not in zero azimuth
c
c------------------------------------------------------------------------
c
c    calls NAG routines...
c
c G05CBF Initialise random number generating routines to give repeatable 
c        sequence
c G05CCF Initialise random number generating routines to give non-repeatable
c        sequence
c G05DDF Pseudo-random real numbers, Normal distribution
c
c E04UDM dummy  CONFUN  routine for use with  E04UCF  when there are no 
c        constraints present
c E04URF Supply optional parameter values to E04UPF
c E04UPF Minimum of a sum of squares, nonlinear constraints, sequential QP
c        method, using function values and optionally 1st derivatives
c
c
c    calls NUMERICAL RECIPES routines...
c
c GAMMQ
c RTBIS
c
c------------------------------------------------------------------------
c
c Inputs:
c   IND1        data files  (unit=4, set in io.inc)
c   IND2        list file   (unit=3, set in io.inc)
c
c Outputs:
c  OUTD1        dcmp output file (unit=7, set in io.inc)
c  OUTD2        impedance output file (unit=8, set in io.inc)
c
c
c------------------------------------------------------------------------
c Notes:
c 
c logical WARM sets whether E04UPF can use Warm Start option. 
c Default is WARM = .FALSE.  This increases run time significantly.
c To speed up the program, change to WARM = .TRUE.
c
c
c------------------------------------------------------------------------
c Country/System Specific Notes:
c
c you must define the FORWARD SLASH and BACK SLASH if it is defined
c differently from ASCII ANSI X3.4
c
c char(47) = UNIX directory separator = BACK SLASH
c char(92) = DOS  directory separator = FORWARD SLASH
c
c========================================================================

      implicit none

c------------------------------------------------------------------------
      real PI
      parameter ( PI = 3.14159265 )
      real MU0
      parameter ( MU0 = 4.*PI*1.E-7 )

      include 'io.inc'
      include 'size.inc'
      include 'ctrlblk.inc'
      include 'version.inc'

c------------------------------------------------------------------------

      complex*16 alpha(4,MAXF,MAXS), a, b, ztmp1(2,2), ztmp2(2,2)

      complex zmeas(2,2,MAXF,MAXS), zrel(2,2,MAXF,MAXS), 
     &        zfit(2,2,MAXF,MAXS), a4, b4

      real*8 alow, ahigh, alowin, ahighin, sazim, sshear, stwist, 
     &       bound, slow, shigh, tlow, thigh, 
     &       e, t, ra, rb, ia, ib,
     &       bu((MAXF*4+2)*MAXS+1), x((MAXF*4+2)*MAXS+1),
     &       bl((MAXF*4+2)*MAXS+1), error1, terr(2,2),
     &       al_dev(4,MAXF,MAXS), 
     &       aa(1,(MAXF*4+2)*MAXS+1), 
     &       r((MAXF*4+2)*MAXS+1,(MAXF*4+2)*MAXS+1),
     &       f(MAXF*8*MAXS), objf, 
     &       fjac(MAXF*8*MAXS,(MAXF*4+2)*MAXS+1),
     &       clamda((MAXF*4+2)*MAXS+1), user(1), cjac(1,1), c(1)

      real std_dev(2,2,MAXF,MAXS), period(MAXF,MAXS), zmax_floor
      real convz2r, convz2p, rd, factor
 
      real permin, permax, per, perband(2, MAXBND),
     &     freqs(MAXF,MAXS), skew(MAXF,MAXS),
     &     rho_a(MAXF,MAXS), rho_b(MAXF,MAXS), 
     &     pha_a(MAXF,MAXS), pha_b(MAXF,MAXS), 
     &     azim_chan(MAXF,MAXS), err(MAXF,MAXS),
     &     shear(MAXF,MAXS), twist(MAXF,MAXS), azim(MAXF,MAXS), 
     &     sk_b(2,MAXF,MAXS), er_b(2,MAXF,MAXS), sh_b(2,MAXF,MAXS), 
     &     tw_b(2,MAXF,MAXS), chiprob, p, chiv, prob, gammq,
     &     az_b(2,MAXF,MAXS), 
     &     rho_a_b(2,MAXF,MAXS), rho_b_b(2,MAXF,MAXS), 
     &     pha_a_b(2,MAXF,MAXS), pha_b_b(2,MAXF,MAXS), 
     &     chan_b(2,MAXF,MAXS),
     &     permin_in, permax_in, azimuth, azimuth1, toterr

      real perdat(MAXDAT)
      real z_floor, rho_floor, pha_floor
      complex zxx(MAXDAT), zxy(MAXDAT), zyx(MAXDAT), zyy(MAXDAT)
      real zxxe(MAXDAT), zxye(MAXDAT), zyxe(MAXDAT), zyye(MAXDAT)

c...these arrays store the final resuults from each band
      real twists(0:MAXREL,MAXFBD,MAXS), 
     &     shears(0:MAXREL,MAXFBD,MAXS),
     &     azims (0:MAXREL,MAXFBD,MAXS),
     &     errs  (0:MAXREL,MAXFBD,MAXS),
     &     skews (0:MAXREL,MAXFBD,MAXS),
     &     rho_as(0:MAXREL,MAXFBD,MAXS), 
     &     rho_bs(0:MAXREL,MAXFBD,MAXS),
     &     pha_as(0:MAXREL,MAXFBD,MAXS), 
     &     pha_bs(0:MAXREL,MAXFBD,MAXS)
      real array(0:MAXREL), var, xmean, tol, bandwidth, bandlap
      real zrat, s, sl, su
      real realz, imagz
      real objstart(MAXBND), objend(MAXBND)

      integer m, i, n, ifail, o, oo, pos, ios, nf, nf1, nf2
      integer k, count(MAXS), nrel, nrel0, jrel, j, nsite, k1, kstart,
     &         j1, j2, jtot, ktot, nsitetot
      integer dpos, spot
      integer spos, tpos, fpos
      integer lda, ldcj, ldfj, ldr, iter, iuser(1)
      integer istate(NVARIABLES), neq
      integer it, it_inc, nu, linsite, nstate
      integer lnblnk
      integer nband, iband, bands(MAXF), isite(MAXS)
      integer iseed, islash
      integer MODE, NROWJ, NEEDC

      integer print_level

      character normtype*6

      character*80 inlist, infile(MAXS),
     &             tsite, string, filein(MAXS), file, line
      character*20 insite(MAXS), insite1
      character*1  yans
      character*11 bdate
    
      logical  STATS               /.FALSE./,
     &         LBOUND              /.TRUE./, 
     &         GRAD                /.FALSE./,
     &         EXPERT              /.FALSE./,
     &         SINGLE              /.FALSE./,
     &         FIRST_TIME          /.TRUE./,
     &         DEBUG,
     &         SCALE,
     &         WARM,
     &         MAX_Z_ERROR_FLOOR   /.TRUE./,
     &         EDI                 /.FALSE./

      character*24 date_string
      character*1 bslash, fslash

      common /cmn1/ alpha, al_dev, normtype
      common /cmn2/nsite, count
      common /cmndbg/ DEBUG

      external objfun, e04udm, gammq, rtbis

      real atand
      real*8 dtand, g05ddf
      
      integer iwork(LW)
      real*8 work(LIW)

      data azim/MAXFS*-9999./
      data istate/NVARIABLES*0/
      
c-------------------------------------------------------------------------------
c...set character set dependent variables

      bslash = char(47)
      fslash = char(92)
      
c-------------------------------------------------------------------------------

c...set build date
      bdate = '04-March-2005'

c...get date
      call fdate(date_string)

c-------------------------------------------------------------------------------


      write(*,*)' '
      write(*,*)'                Welcome to strike'
      write(*,*)'     Written by Alan G Jones & Gary W McNeice'
      write(*,*)' '
      write(*,*)'                  version: ', ver,'b'
      write(*,*)'       *****THIS IS STILL A BETA TEST VERSION*****'
      write(*,*)'              Build date: ', bdate
      write(*,*)'        Todays date: ', date_string
      write(*,*)' '
      write(*,*)'  Version 4.+ reads in edi or J-format files'
      write(*,*)' '
      write(*,*)'Note: field units assumed for EDI file, not S.I.'
      write(*,*)' '
      write(*,*)'               Version limits:'
      write(*,*)'      Max. no. sites:       ', MAXS
      write(*,*)'      Max. no. freqs:       ', MAXF
      write(*,*)'      Max. no. bands :      ', MAXBND
      write(*,*)'      Max. no. freqs/band : ', MAXFBD
      write(*,*)'      Max. no. data :       ', MAXDAT
      write(*,*)'      Max. no. realizations:', MAXREL
      write(*,*)' '
      write(*,*)' Min. no. realizations for statistics: ', MINREL
      write(*,*)' '


      nrel0 = 0
      permin_in = 1.e-4
      permax_in = 1.e+4
      print_level = 0

c...set Optimality tolerance
      tol = 0.0001

c...set starting point
      sazim  = 20.0D+00
      sshear = 10.0D+00
      stwist = 20.0D+00

      rd = 180.0/PI

c...request if EXPERT
      yans = 'N'
      call yin( 'Expert mode (y/n)', yans )
      if( yans.eq.'Y' ) then
        EXPERT = .TRUE.
      else
        EXPERT = .FALSE.
      endif


      if( EXPERT ) then
      
c...ask for general print level
      iprint = 1
      call iin( 'Give print level (0/1/2/3/10/20/30)', iprint )

c...ask for DEBUG mode
        yans = 'Y'
        call yin( 'Debug mode (y/n)', yans )
        if( yans.eq.'Y' ) then
          DEBUG = .TRUE.
          write(*,*)'Debug mode writes out the estimates from the',
     &              ' to various files:'
          open(unit=30,file='azimuths.str')
          open(unit=31,file='rho_as.str')
          open(unit=32,file='rho_bs.str')
          open(unit=33,file='pha_as.str')
          open(unit=34,file='pha_bs.str')
          open(unit=35,file='shears.str')
          open(unit=36,file='twists.str')
          open(unit=37,file='errors.str')
          open(unit=38,file='rels.str')
          open(unit=39,file='impedances.str')
          write(*,*)'impeds.str  : impedances input'
          write(*,*)'azimuths.str: azimuth     estimates'
          write(*,*)'rho_as.str  : rho_a       estimates'
          write(*,*)'rho_bs.str  : rho_b       estimates'
          write(*,*)'pha_as.str  : pha_a       estimates'
          write(*,*)'pha_bs.str  : pha_b       estimates'
          write(*,*)'shears.str  : shear       estimates'
          write(*,*)'twists.str  : twist       estimates'
          write(*,*)'rels.str    : realization estimates'

        else
          DEBUG = .FALSE.
          iprint = 0
        endif

c...ask for print level from NAG routines
        call iin('NAG print level (0, 1, 5, 10, 20 or 30)',print_level)
        if( print_level.gt.0 ) then
          call e04urf('List')
        else
          call e04urf('Nolist')
        endif
        
c...ask if gradients are to be verified
        call lin('Verify gradients (T or F)', GRAD )
        if( GRAD ) then
          call e04urf('Verify level = 3')
        end if

c...ask for optimality tolerance
        call rin( 'Give Optimality tolerance', tol )


c...ask for starting solution
        call rin8('Give starting azimuth', sazim )
        call rin8('Give starting shear', sshear )
        call rin8('Give starting twist', stwist )

c...request whether to fit to scaled impedances
        yans = 'N'
        call yin( 'Fit scaled impedances (Z/sqrt(omega*mu))', yans )
        if( yans.eq.'Y' ) then
          SCALE = .TRUE.
        else
          SCALE = .FALSE.
        endif

c...request whether to use Warm Start option on successive calls to E04UPF
        yans = 'N'
        call yin( 'Warm start successive calls to E04UPF', yans )
        if( yans.eq.'Y' ) then
          WARM = .TRUE.
        else
          WARM = .FALSE.
        endif

c...request whether error floor is based on maximum value
        yans = 'Y'
        call yin( 'Base error floor on maximum impedance', yans )
        if( yans.eq.'Y' ) then
          MAX_Z_ERROR_FLOOR = .TRUE.
        else
          MAX_Z_ERROR_FLOOR = .FALSE.
        endif


        write(*,*)' '

      else

        call e04urf('Nolist')
        SCALE = .FALSE.
        WARM  = .FALSE.

      endif
      
      
c...set NAG optimality tolerance
      string = 'Optimality tolerance = '
      write(string(24:34),'(f10.8)') tol
      call e04urf( string )

c...set NAG print level
      string = 'Major print level = '
      write(string(21:22),'(i2)') print_level
      call e04urf( string )



c-------------------------------------------------------------------------------
c--- Find sites to operate on

      inlist = 'site.lis'
      call tin('Site list (.dat or .edi for single site)?', inlist)
      if( print_level.ge.1 ) then
        write(*,*)'Read in site list: ', inlist(1:lnblnk(inlist))
      endif

      if( index(inlist,'.dat').eq.0 .and. 
     &    index(inlist,'.edi').eq.0 )then
        SINGLE = .FALSE.
        open (unit=IND2, file=inlist, status='old', err=9998)

        do o = 1, MAXS+1
          read(IND2,'(a)', err=100, end=100) line
	  if( print_level.ge.1 ) then
            write(*,*)'Read in filename: ', line(1:lnblnk(line))
	  endif

          if( o.gt.MAXS ) then
            write(*,*)'Too many sites for this version....'
            write(*,*)'nsitetot =', nsitetot,'    MAXS =', MAXS
            stop
          endif

          infile(o) = line

          write(*,'(a5,i2,a9,a80)')'Site ', o,'  file = ',
     &                           infile(o)(1:lnblnk(infile(o)))
        enddo

100     nsitetot = o-1
        close(unit=17)

      else
        SINGLE = .TRUE.
        infile(1) = inlist
        nsite = 1
        nsitetot = 1
      endif
      write(*,'(a,i2)')'No. of sites to process =', nsitetot

c...read in azimuth information from first (only) file
      o = 1
      write(*,*)'Opening first file =',
     &                         infile(1)(1:lnblnk(infile(1)))
      open(unit=IND1, file=infile(1), status='old', err=9999)
      
c...check if EDI file
      if( index(infile(1),'.edi').ne.0 ) then
	EDI = .TRUE.
      else
        EDI = .FALSE.
      endif
      
      if( EDI ) then
	call ediazimuth( ind1, azimuth1 )
      else
        line = ' '
        do while(line(1:8).ne.'>AZIMUTH')
          read(IND1,'(a)',end=101) line
        enddo
        read(line(13:),*) azimuth1
        goto 102
101     write(*,*)'AZIMUTH not found in J-format file >',
     &             infile(1)(1:lnblnk(infile(1)))
        write(*,*)'...set to zero'
        azimuth1 = 0.0
      endif
102   write(*,*)'Azimuth from first file =', azimuth1
      close(unit=IND1)


c...request for error floor level
      z_floor = 1.75
      if( MAX_Z_ERROR_FLOOR ) then
        call rin( 'Give impedance relative error floor of maximum'//
     &            ' value (in %)', z_floor )
      else
        call rin( 'Give impedance relative error floor of individual'//
     &            ' value (in %)', z_floor )
      endif
c...convert z_floor from percentage to fractional floor
      z_floor = z_floor/100.
      pha_floor = atand( z_floor )
      rho_floor = 100.*(2.*z_floor + z_floor**2.)
      write(*,'(a,f5.2)')'Error floor: z_floor (%)   =',z_floor*100.
      write(*,'(a,f5.2)')'equivalent pha_floor (deg) =',pha_floor
      write(*,'(a,f5.2)')'equivalent rho_floor (%)   =',rho_floor



c...request normalization type
103   normtype = 'GAVSD2'
      call tin( 'Give normalization type (? for list)', normtype )
104   if( normtype.eq.'?' ) then
        write(*,*)'L2          L2 unweighted normalization'
        write(*,*)'MAXSD       maximum s.d. weighted normalization'
        write(*,*)'GAVSD       geom. av. s.d. weighted normalization'
        write(*,*)'GAVSD2      geom. av. s.d. (Zxx,Zyy) & (Zxy,Zyx) '//
     &                         'weighted normalization'
        write(*,*)'SUMSQ       sqrt(sum squares s.d.) '//
     &                         'weighted normalization'
        write(*,*)'SUMSQ2      sqrt(sum squares s.d. (Zxx,Zyy) & '//
     &                         '(Zxy,Zyx)) weighted normalization'
        goto 103
      endif
      if( normtype(1:2).ne.'L2'     .and.
     &    normtype(1:5).ne.'MAXSD'  .and.
     &    normtype(1:5).ne.'GAVSD'  .and.
     &    normtype(1:6).ne.'GAVSD2' .and.
     &    normtype(1:6).ne.'SUMSQ'  .and.
     &    normtype(1:6).ne.'SUMSQ2' ) then
        write(*,*)'Unacceptable response: give one of the following'
        normtype = '?'
        goto 104
      endif
      if( DEBUG ) then
        write(*,*)'normtype entered = ', normtype
      endif



c...request frequency bounds
      call rin('Enter minimum period?', permin_in )
      call rin('Enter maximum period?', permax_in )
      bandwidth = log10(permax_in/permin_in)
      call rin('Enter bandwidth (no. of period decades)?', bandwidth)
      bandlap = 0.0
      call rin('Enter overlap   (no. of period decades)?', bandlap)
      nband = nint( log10(permax_in/permin_in)/(bandwidth-bandlap) )
      if( nband.eq.0 ) nband = 1
      write(*,'(a,i2)')'No. of bands to analyze =', nband

c...check if number of bands too large
      if( nband.gt.MAXBND ) then
        write(*,*)' '
        write(*,*)'WARNING: NUMBER OF BANDS TOO LARGE FOR THIS VERSION'
        write(*,*)' '
        write(*,*)'nband reset to MAXBND =', MAXBND
        nband = MAXBND
      endif

      do iband = 1, nband
        permin = 10.**(log10(permin_in) + float(iband-1)
     &                                       *(bandwidth-bandlap))
        permax = 0.999*10.**(log10(permin) + bandwidth)
        if( permin.lt.1. .and. permax.le.1. ) then
          write(*,'(a,i2,a,f10.2,a,f10.2,a)')'Band ', iband, 
     &          ':  Periods ', 1./permin,' Hz    -  ', 1./permax,' Hz'
        else if( permin.lt.1. .and. permax.gt.1. ) then
          write(*,'(a,i2,a,f10.2,a,f10.2,a)')'Band ', iband, 
     &          ':  Periods ', 1./permin,' Hz    -  ', permax,' s'
        else 
          write(*,'(a,i2,a,f10.2,a,f10.2,a)')'Band ', iband, 
     &          ':  Periods ', permin,' s     -  ', permax,' s'
        endif
      enddo


c...ask if bounds are needed  A PROBLEM WITH THIS - TO BE FIXED LATER
      if( EXPERT ) then
        yans = 'Y'
c        call yin('Place bounds on parameters (Y or N) ?', yans )
        if( yans.eq.'Y' ) then
          LBOUND = .TRUE.
        else
          LBOUND = .FALSE.
        endif
      endif

c...get bounds
      if( LBOUND ) then
        alowin  = -1080.0D+00
        ahighin =  1080.0D+00
        slow  = -45.0D+00
        shigh =  45.0D+00
        tlow  = -60.0D+00
        thigh =  60.0D+00

        write(*,*) ' '
        write(*,*) ' Regional azimuth bounds ', sngl(alowin),' ',
     &                                          sngl(ahighin)
        write(*,*) '            Shear bounds ', sngl(slow),' ',
     &                                          sngl(shigh)
        write(*,*) '            Twist bounds ', sngl(tlow),' ',
     &                                          sngl(thigh)
        write(*,*) ' '

        yans = 'N'
        call yin('Change bounds from standard bounds (y/n)?',yans)
        if( yans.eq.'Y' ) then
          write(*,*) ' '
          call rin8('Give strike lower bound', alowin )
          call rin8('Give strike upper bound', ahighin )
          call rin8('Give shear lower bound',  slow )
          call rin8('Give shear upper bound',  shigh )
          call rin8('Give twist lower bound',  tlow )
          call rin8('Give twist upper bound',  thigh )
          write(*,*) ' '
        endif
      endif

c...are statistics needed?
      yans = 'N'
      call yin( 'Do statistics (y/n)', yans )
      if( yans.eq.'Y' ) then
        STATS = .TRUE.
      else
        STATS = .FALSE.
      endif

      if( STATS ) then
        nrel = MAXREL
        call iin( 'Give number of realizations (default is maximum'//
     &            ' permitted)', nrel )
        if( EXPERT ) then
          yans = 'Y'
          call yin('Random starting seed (y/n)?', yans )
          if( yans.eq.'Y' ) then
            call g05ccf
          else
            iseed = 23427263
            call iin( 'Give starting seed integer', iseed )
            call g05cbf(iseed)
          endif
        else
          call g05ccf
        endif
      else
        nrel = 0
      endif

c===============================================================================
c      BIG LOOP OVER nband BANDS
c
c      This opens each data file and reads in the relevant data for
c      the period band of interest. Then the datafile is closed.

      kstart = 0
      do 1000 iband = 1, nband

        permin = 10.**(log10(permin_in) + float(iband-1)
     &                                       *(bandwidth-bandlap))
        permax = 0.999*10.**(log10(permin) + bandwidth)

        if( DEBUG ) then
          write(39,*)'band, permin, permax =', iband, permin, permax
        endif

        perband(1,iband) = permin
        perband(2,iband) = permax

        write(*,*)' '
        write(*,'(a,i2,a,f12.4,a,f12.4)') 'Doing band ', iband, 
     &             '    Period range',
     &             permin, '  -', permax

c...open each file and check that data exist for the band
c   find start and end freq to read in

        if( .not.SINGLE ) then
          close(unit=IND2)
          if( DEBUG ) then
            write(*,*)'Opening file: ', inlist(1:lnblnk(inlist))
          endif
          open (unit=IND2, file=inlist, status='old', err=9999)
        endif

        nsite = 0
        do 110 o = 1, nsitetot
          nf = 0
          if( .not.SINGLE ) then
            read(IND2,'(a)',err=666) infile(o)
          endif

c...check if EDI file
          if( index(infile(o),'.edi').ne.0 ) then
	    EDI = .TRUE.
	  else
	    EDI = .FALSE.
	  endif	 

c...pick off site name
          if( EDI ) then	  
            fpos = index(infile(o),'.edi') + 3
	  else
            fpos = index(infile(o),'.dat') + 3
	  endif
          tpos = index(infile(o),'/')
          if( tpos .eq. 0 ) then
            spos = 1
          else
            spos = fpos - 1
            do j = 1, fpos-3
              spos = spos - 1
              if(infile(o)(spos:spos) .eq. '/') then
                spos = spos + 1
                goto 194
              endif
            enddo
          endif
194       insite1 = infile(o)(spos:fpos)


c...open file
          if( DEBUG ) then
            write(*,*)'Opening file: ',infile(o)
          endif
          close(unit=IND1)
          open(unit=IND1, file=infile(o), status='old', err=9999)

c...get number of periods in file
          if( EDI ) then
	    call ediformatfin( IND1, permin, permax, nf )
	  else
	    call   jformatfin( IND1, permin, permax, nf )
	  endif
	  
          write(*,'(3a,i3)')'No. freqs for site...',
     &                          insite1, ':', nf

          if( DEBUG ) then
            write(39,*)'Site =', o, 
     &                 '  file =', infile(o)(1:lnblnk(infile(o))),
     &                 '  no freqs =', nf
          endif

          if( nf.gt.MAXFBD ) then
            write(*,*)'Too many frequencies in this band for this '//
     &                'version'
            write(*,*)'MAXFBD = ', MAXFBD
            write(*,*)'Increase this value in your size.inc file'
            stop
          endif

          if( nf.gt.0 ) then
            nsite = nsite+1
            filein(nsite) = infile(o)
            isite(nsite) = o
          endif

110     close(unit=IND1)

        write(*,'(a,i2)')'No. of sites for this band....', nsite
        if( nsite.eq.0 ) then
          write(*,'(a)')'No data to analyze - skipping to next band'
          goto 1000
        endif

c...open site          

      do o = 1, nsite
        if( SINGLE ) then
          file = inlist
        else
          file = filein(o)
        endif
        if( EDI ) then
          fpos = index(file,'.edi') + 3
        else
          fpos = index(file,'.dat') + 3
        endif
        tpos = index(file,'/')
        tsite = file
        if ( tpos .eq. 0 ) then
          spos = 1
        else
          spos = fpos - 1
          do j = 1, 80
            spos = spos - 1
            if (tsite(spos:spos) .eq. '/') then
              spos = spos + 1
              goto 195
            endif
          enddo
        endif
195     if( DEBUG ) then
          write(*,'(3a)')'Reading file > ', tsite(spos:fpos),'...'
        endif


c...open file
        close(unit=IND1)
        open(unit=IND1, file=file, status='old', err=9999)

c...read in azimuth information - check that all sites have same azimuth
        if( EDI ) then
	  call ediazimuth( ind1, azimuth )
	else
          line = ' '
          do while(line(1:8).ne.'>AZIMUTH')
            read(IND1,'(a)',end=200) line
          enddo
          read(line(13:),*) azimuth
	endif
        goto 205

200     write(*,'(a)')'No azimuth information in file: set to zero'
        azimuth = 0.

205     if( azimuth.ne.azimuth1 ) then
          write(*,*)'ERROR: files have different azimuths:'
          write(*,*)'Previous azimuth:', azimuth1
          write(*,*)'This     azimuth:', azimuth
	  write(*,*)'Sorry - not permitted in this version of strike'
          stop
        endif


210     rewind(unit=IND1)
        alow  = alowin  - azimuth
        ahigh = ahighin - azimuth

c...pick off site name
        insite(o) = tsite(spos:fpos)

c...read in impedances and impedance estimates from J-format or EDI-format file

        if( EDI ) then
          call edizdatin( IND1, n, perdat, zxx, zxxe, zxy, zxye, 
     &                                     zyx, zyxe, zyy, zyye )
          if( DEBUG ) write(*,*)'strike: returned from edizdatin'
        else
          call   jzdatin( IND1, n, perdat, zxx, zxxe, zxy, zxye, 
     &                                     zyx, zyxe, zyy, zyye )
          if( DEBUG ) write(*,*)'strike: returned from jzdatin'
        endif
        if( DEBUG ) write(*,*)'strike: Read in data - n =', n

        count(o) = 0
        nf1 = 0
        nf2 = 0
        do 300 m = 1, n
          per = perdat(m)

          if( per.ge.permin .and. nf1.eq.0 ) nf1 = m
          if( per.gt.permax .and. nf2.eq.0 ) nf2 = m-1

          if( (per.lt.permin) .or. (per.gt.permax) ) goto 300

          count(o) = count(o) + 1
          if( count(o).gt.MAXF ) then
            write(*,*)' '
            write(*,*)'***** Too many frequencies for this version'
            write(*,*)'MAXF =', MAXF,'   count(o) =', count(o)
            write(*,*)'     only the first ',MAXF,' will be used'
            write(*,*)' '
            count(o) = MAXF
            goto 310
          endif

          if( SCALE ) then
            factor = sqrt(per/(MU0*2.*PI))
          else
            factor = 1.
          endif

          period(count(o),o)      = per
          zmeas(1,1,count(o),o)   = factor*zxx(m)
          zmeas(1,2,count(o),o)   = factor*zxy(m)
          zmeas(2,1,count(o),o)   = factor*zyx(m)
          zmeas(2,2,count(o),o)   = factor*zyy(m)
          std_dev(1,1,count(o),o) = factor*zxxe(m)
          std_dev(1,2,count(o),o) = factor*zxye(m)
          std_dev(2,1,count(o),o) = factor*zyxe(m)
          std_dev(2,2,count(o),o) = factor*zyye(m)

          if( DEBUG ) then
            write(39,*)'Period =', per
            write(39,*) zmeas(1,1,count(o),o), zmeas(1,2,count(o),o)
            write(39,*) zmeas(2,1,count(o),o), zmeas(2,2,count(o),o)
          endif

c...check error > error floor. If smaller, set to error floor
          if( MAX_Z_ERROR_FLOOR ) then
            zmax_floor = z_floor*max( 
     &                 cabs(zmeas(1,1,count(o),o)),
     &                 cabs(zmeas(1,2,count(o),o)),
     &                 cabs(zmeas(2,1,count(o),o)),
     &                 cabs(zmeas(2,2,count(o),o)) )
            do i = 1, 2
              do j = 1, 2
                std_dev(i,j,count(o),o) = max( 
     &               std_dev(i,j,count(o),o), zmax_floor )
              enddo
            enddo
          else
            do i = 1, 2
              do j = 1, 2
                std_dev(i,j,count(o),o) =  
     &                    max(std_dev(i,j,count(o),o),
     &              z_floor*cabs(zmeas(i,j,count(o),o)) )
              enddo
            enddo
          endif

300     continue

310     if( o.eq.1 ) then
          bands(iband) = count(o)
          write(*,'(a,i4)') 'No. of frequencies in this band =', 
     &                           count(o)
          if( count(o).eq.0 ) goto 1000
        endif
      enddo
666   if( SINGLE ) then
        close(4)
      else
        close(17)
      endif

c...all data read in: no. of sites = nsite
      nsite = o - 1
      if( DEBUG ) then
        write(*,*) '... ',nsite, ' sites read in ...'
      endif

      if( nrel.gt.MAXREL ) then
        nrel = MAXREL
        write(*,*) ' '
        write(*,*) 'no. of realizations reset to', MAXREL
      endif


c...calculate the number of unknowns (parameters) - n

c********* real and imag. imp.   +  sh & tw + strike
      n = 0
      neq = 0
      do o = 1, nsite
        n = n + count(o)*4  + 2
        neq = neq + count(o)*8
      end do
      n = n + 1
      lda = 1
      ldcj = 1
      ldfj = (MAXF*8)*MAXS
      ldr = n


c...calculate alpha standard deviation
       do o = 1, nsite
        do k = 1, count(o)
c...v3.0
c         al_dev(1,k,o) = sqrt(std_dev(1,1,k,o)**2+std_dev(2,2,k,o)**2)
c         al_dev(2,k,o) = sqrt(std_dev(1,2,k,o)**2+std_dev(2,1,k,o)**2)
c...v3.1
          al_dev(1,k,o) = std_dev(1,1,k,o)
          al_dev(2,k,o) = std_dev(1,2,k,o)
          al_dev(3,k,o) = std_dev(2,1,k,o)
          al_dev(4,k,o) = std_dev(2,2,k,o)
        enddo
       enddo

c...............................................................................
c...BIG LOOP OVER nrel REALIZATIONS.............................................

      do 500 jrel = 0, nrel
        write(*,'(a,i2)')'Doing realization...', jrel
        if( DEBUG ) then
          write(38,'(a,i2)')'Realization...', jrel
        endif
        do o = 1, nsite
          do k = 1, count(o)
            if( DEBUG ) write(38,*)'site =',o,'  freq=',k
            do j = 1, 2
              do i = 1, 2
                if( jrel.eq.0 ) then
                  realz = 0.
                  imagz = 0.
                else
                  realz=sngl(g05ddf(0.d0,dble(std_dev(i,j,k,o))))
                  imagz=sngl(g05ddf(0.d0,dble(std_dev(i,j,k,o))))
                endif
                zrel(i,j,k,o) = zmeas(i,j,k,o) + cmplx(realz, imagz)
                if( DEBUG ) then
                  write(38,*)i,j, 
     &                      real(zmeas(i,j,k,o)),aimag(zmeas(i,j,k,o)),
     &                      std_dev(i,j,k,o), realz, imagz
                 endif
              enddo
            enddo
          enddo
        enddo

c...calculate alphas for all frequencies and sites
c   Eqns 33 of GB89
      do o = 1, nsite
        do k = 1, count(o)
          alpha(1,k,o) = (zrel(1,1,k,o) + zrel(2,2,k,o))
          alpha(2,k,o) = (zrel(1,2,k,o) + zrel(2,1,k,o))
          alpha(3,k,o) = (zrel(2,1,k,o) - zrel(1,2,k,o))           
          alpha(4,k,o) = (zrel(1,1,k,o) - zrel(2,2,k,o))
        end do
      end do

c...set up bounds

      if( .not.LBOUND ) goto 777

c---( (no.freqs.*4 + 2)  * no.sites + 1 ) unknowns:
c********* (real and imag. imp.+sh & tw)*no.sites + strike
c   1: Re(a) freq. 1 site 1
c   2: Im(a) freq. 1 site 1
c   3: Re(b) freq. 1 site 1
c   4: Im(b) freq. 1 site 1
c   5: Re(a) freq. 2 site 1
c   6: Im(a) freq. 2 site 1
c   7: Re(b) freq. 2 site 1
c   8: Im(b) freq. 2 site 1
c   .
c   .
c   .
c   (no. freqs. * 4) + 1: t for site 1
c   (no. freqs. * 4) + 2: e for site 1
c   same for other sites
c   .
c   .
c   .
c   azimuth for all sites

c...set impedance bound

      bound = 0.0D+00
      do o = 1, nsite
       do k = 1, count(o)
        if(cabs(zrel(1,1,k,o)).gt.bound)bound=cabs(zrel(1,1,k,o))
        if(cabs(zrel(1,2,k,o)).gt.bound)bound=cabs(zrel(1,2,k,o))
        if(cabs(zrel(2,1,k,o)).gt.bound)bound=cabs(zrel(2,1,k,o))
        if(cabs(zrel(2,2,k,o)).gt.bound)bound=cabs(zrel(2,2,k,o))
       enddo
      enddo
c...set impedance bound to 4x largest value found
      bound = bound * 4.0D0


c...each site

      pos = 0
      do o = 1, nsite
        if (o.gt.1) then
          pos = pos + count(o-1)*4+2
        end if
        do k = 1, count(o)
          bl(pos+(k*4)-3) = -bound
          bu(pos+(k*4)-3) =  bound
          bl(pos+(k*4)-2) = -bound
          bu(pos+(k*4)-2) =  bound
          bl(pos+(k*4)-1) = -bound
          bu(pos+(k*4)-1) =  bound
          bl(pos+k*4)     = -bound
          bu(pos+k*4)     =  bound
        enddo
        bl(pos+(count(o)*4)+1) = dtand(tlow)
        bu(pos+(count(o)*4)+1) = dtand(thigh)
        bl(pos+(count(o)*4)+2) = dtand(slow)
        bu(pos+(count(o)*4)+2) = dtand(shigh)
      enddo

c...set strike bounds
      bl(n) = alow  / rd  
      bu(n) = ahigh / rd  
 
c...set starting point
c...each site

777   pos = 0
      do o = 1, nsite
        if (o.gt.1) then
          pos = pos + count(o-1)*4+2
        endif
        do k = 1, count(o)
          x(pos+(k*4)-3) = real(zrel(1,2,k,o))
          x(pos+(k*4)-2) = imag(zrel(1,2,k,o))
          x(pos+(k*4)-1) = 0.8d0 * real(zrel(1,2,k,o))
          x(pos+k*4)     = 0.8d0 * imag(zrel(1,2,k,o))
        enddo
        x(pos+(count(o)*4)+1) = dtand(stwist)
        x(pos+(count(o)*4)+2) = dtand(sshear)
      enddo
c...strike
      x(n) = sazim / rd

c...call function to get starting value
      nstate = 1
      if( DEBUG .and. iprint.ge.5 ) then
        write(*,*)'Calling objfun first time with...'
        write(*,*)'neq    =', neq
        write(*,*)'n      =', n
        write(*,*)'ldfj   =', ldfj
        write(*,*)'nstate =', nstate
        write(*,*)'x(1:n) =>'
        do i = 1, n
          write(*,*) i, x(i)
        enddo
      endif
      call objfun(2,neq,n,ldfj,x,f,fjac,nstate,iuser,user)
      objf = 0.D0
      do i = 1, neq
       objf = objf + f(i)**2.d0
      enddo
      objstart(iband) = sngl(objf)
      if( DEBUG ) then
        write(*,'(a,g12.4)')'Objective function total =',
     &                                     sngl(objf)
        if( iprint.ge.5 ) then
          do i = 1, neq
            write(*,*)'Equation', i, '   f =', f(i)
          enddo
        endif
      endif


      it = max(50,25*n)
      it_inc = 0
      if( EXPERT .and. jrel.eq.0 .and. iband.eq.1 ) then
        call iin('Maximum number of iterations', it)
      endif

c...set ifail=-1
800   ifail = -1
      iter = 0
      if( .not.WARM ) then
        do i = 1, n
          istate(i) = 0
        enddo
      endif

      if( DEBUG .and. iprint.ge.5 ) then
        write(*,*)'Entering E04UPF with...'
          write(*,*) ' 1: m (neq)      = ', neq
          write(*,*) ' 2: n (n)        = ', n
          write(*,*) ' 3: NCLIN        = ', NCLIN
          write(*,*) ' 4: NCNLN        = ', NCNLN
          write(*,*) ' 5: lda          = ', lda
          write(*,*) ' 6: ldcj         = ', ldcj
          write(*,*) ' 7: ldfj         = ', ldfj
          write(*,*) ' 8: ldr          = ', ldr
          write(*,*) ' 9: a (aa)       = NOT REFERENCED'
          write(*,*) '10: bl           = lower bounds'
          write(*,*) '11: bu           = upper bounds'
          write(*,*) '12: CONFUN       = E04UDM'
          write(*,*) '13: OBJFUN       = objfun'
          write(*,*) '14: iter         = ', iter
          write(*,*) '15: istate(1)    =', istate(1)
          write(*,*) '    WARM         =', WARM
          write(*,*) '16: c            = output array'
          write(*,*) '17: cjac         = output array'
          write(*,*) '18: f            = output array'
          write(*,*) '19: fjac         = output array'
          write(*,*) '20: clambda      = output array'
          write(*,*) '21: objf         = output parameter'
          write(*,*) '22: r            = output array'
          write(*,*) '23: x            = initial solution'
          write(*,*) '24: iwork        = integer workspace array'
          write(*,*) '25: liwork (LIW) = ', LIW
          write(*,*) '26: work         = real workspace array'
          write(*,*) '27: lwork (LW)   = ', LW
          write(*,*) '28: iuser        = iuser(1)'
          write(*,*) '29: user         = user(1)'
          write(*,*) '30: ifail        = ', ifail

          do i = 1, n
            write(*,*) i, sngl(bl(i)), sngl(x(i)), sngl(bu(i))
          enddo
      endif


      call e04udm(MODE,NCNLN,N,NROWJ,NEEDC,X,C,CJAC,NSTATE,IUSER,
     *                  USER)
      call e04upf(neq,n,NCLIN,NCNLN,lda,ldcj,ldfj,ldr,aa,bl,bu,
     &             e04udm,objfun,iter,istate,c,cjac,f,fjac,clamda,objf,
     &             r,x,iwork,LIW,work,LW,iuser,user,ifail)
      objend(iband) = sngl(objf)

      if( DEBUG ) write(*,*)'Returning from E04UPF...ifail =',ifail

c*******************************************************
c
c   test ifail

      if( ifail.ne.0 ) then
        write(*,*)'error return from e04upf: ifail =', ifail

c...ifail=1: accuracy achieved, but sequence not converged: treat as soft error
        if( ifail.eq.1 ) then
          if( DEBUG ) then
            write(*,*) 'ifail = 1, check one of the following'
            write(*,*)'n      = ', n
            write(*,*)'bl     = ', (bl(i),i=1,n)
            write(*,*)'bu     = ', (bu(i),i=1,n)
            write(*,*)'liwork    = ', liw
            write(*,*)'lwork     = ', lw
          endif

c...ifail=4: no. of iterations been reached. 
c           Increase no. of iterations twice
        elseif( ifail.eq.4 ) then
          it_inc = it_inc + 1
          if( it_inc.le.2 ) then 
            it = 5*it
            write(*,*)'ifail = 4; increasing iterations =', it
            goto 800
          endif


c...ifail=6: does not satisfy first-order Kuhn-Tucker conditions
        elseif( ifail.eq.6 ) then
            write(*,*)'ifile = 6; soln. possibly suspect'


c...ifail=9: invalid input parameter
        elseif( ifail.eq.9 ) then
          write(*,*) 'ifail = 9, invalid input parameter'
          write(*,*) ' 1: m (neq)      = ', neq
          write(*,*) ' 2: n (n)        = ', n
          write(*,*) ' 3: NCLIN        = ', NCLIN
          write(*,*) ' 4: NCNLN        = ', NCNLN
          write(*,*) ' 5: lda          = ', lda
          write(*,*) ' 6: ldcj         = ', ldcj
          write(*,*) ' 7: ldfj         = ', ldfj
          write(*,*) ' 8: ldr          = ', ldr
          write(*,*) ' 9: a (aa)       = NOT REFERENCED'
          write(*,*) '10: bl           = lower bounds'
          write(*,*) '11: bu           = upper bounds'
          write(*,*) '12: CONFUN       = E04UDM'
          write(*,*) '13: OBJFUN       = objfun'
          write(*,*) '14: iter         = ', iter
          write(*,*) '15: istate(1)    =', istate(1)
          write(*,*) '    WARM         =', WARM
          write(*,*) '16: c            = output array'
          write(*,*) '17: cjac         = output array'
          write(*,*) '18: f            = output array'
          write(*,*) '19: fjac         = output array'
          write(*,*) '20: clambda      = output array'
          write(*,*) '21: objf         = output parameter', objf
          write(*,*) '22: r            = output array'
          write(*,*) '23: x            = initial solution'
          write(*,*) '24: iwork        = integer workspace array'
          write(*,*) '25: liwork (LIW) = ', LIW
          write(*,*) '26: work         = real workspace array'
          write(*,*) '27: lwork (LW)   = ', LW
          write(*,*) '28: iuser        = iuser(1)'
          write(*,*) '29: user         = user(1)'
          write(*,*) '30: ifail        = ', ifail

          do i = 1, n
            write(*,*) i, sngl(bl(i)), sngl(x(i)), sngl(bu(i))
          enddo
          stop


c...other ifails are fatal
        else
          write(*,*)'IFAIL condition treated as fatal: contact Alan'
          stop
        endif

      endif

      if( DEBUG .and. jrel.eq.0 ) then
        write(*,*) ' '
        write(*,*)'ifail after e04upf = ', ifail
        write(*,*)'    function value = ', sngl(objf)
        write(*,*) ' '
      endif
 
c...loop over sites                   

      spot = 0
      do 400 o = 1, nsite
        oo = isite(o)
        if( EDI ) then
          dpos = index(filein(o),'.edi')
        else
          dpos = index(filein(o),'.dat')
        endif
        tsite = insite(o)

c...calculate estimated impedances

        if( o.gt.1 ) then
          spot = spot + count(o-1)*4 + 2
        endif
        t = x(spot + count(o)*4 + 1)
        e = x(spot + count(o)*4 + 2)

c...loop over frequencies
        do k = 1, count(o)

c.....calculate "local" freq count and check against bounds
          k1 = kstart + k
          if( k1.gt.MAXFBD ) then
            write(*,*)'EXCEEDING BOUNDS FOR NUMBER OF FREQUENCIES'
            write(*,*)'Re-compile with MAXFBD larger'
            stop
          endif

          if( SCALE ) then
            factor = sqrt(period(k,o)/(MU0*2.*PI))
          else
            factor = 1.
          endif


          ra = x(spot+k*4-3)
          ia = x(spot+k*4-2)
          rb = x(spot+k*4-1)
          ib = x(spot+k*4)
c         write(*,*) k,ra,ia,rb,ib
          a  = dcmplx(ra,ia)
          b  = dcmplx(rb,ib)
          a4 = cmplx(sngl(ra),sngl(ia))/factor
          b4 = cmplx(sngl(rb),sngl(ib))/factor
          call estim_imp(ztmp2,e,t,a,b,x(n))
          zfit(1,1,k,o) = ztmp2(1,1)
          zfit(1,2,k,o) = ztmp2(1,2)
          zfit(2,1,k,o) = ztmp2(2,1)
          zfit(2,2,k,o) = ztmp2(2,2)

c...calculate error of impedance fit
          ztmp1(1,1) = zmeas(1,1,k,o)
          ztmp1(1,2) = zmeas(1,2,k,o)
          ztmp1(2,1) = zmeas(2,1,k,o)
          ztmp1(2,2) = zmeas(2,2,k,o)
          terr(1,1)  = std_dev(1,1,k,o)
          terr(1,2)  = std_dev(1,2,k,o)
          terr(2,1)  = std_dev(2,1,k,o)
          terr(2,2)  = std_dev(2,2,k,o)
          call calc_error(ztmp1,ztmp2,terr, error1)
          errs(jrel,k1,oo) = sngl(error1)

c...calculate apparent resistivities and phases

          rho_as(jrel,k1,oo) = convz2r(a4,period(k,o))
          rho_bs(jrel,k1,oo) = convz2r(b4,period(k,o))
          pha_as(jrel,k1,oo) = convz2p(a4,period(k,o))
          pha_bs(jrel,k1,oo) = convz2p(b4,period(k,o))

c...calculate frequency
          freqs(k1,oo) = period(k,o)

c...calculate azimuth, twist, shear
          azim (k,oo) = sngl(x(n)*rd)   + azimuth
          shear(k,oo) = atand(sngl(e))
          twist(k,oo) = atand(sngl(t))

c...set azimuth in -180 to +180 range
          do while( azim(k,oo).lt.-180. )
            azim(k,oo) = azim(k,oo) + 360.
          enddo
          do while( azim(k,oo).gt.180. )
            azim(k,oo) = azim(k,oo) - 360.
          enddo

          azims (jrel,k1,oo) = azim (k,oo)
          shears(jrel,k1,oo) = shear(k,oo)
          twists(jrel,k1,oo) = twist(k,oo)

c...calculate channelling azimuth
          if(abs(shear(k,oo)).gt.20.0)then
            if(shear(k,oo).gt.0.0)then
              azim_chan(k1,oo)=45.0+azim(k,oo)+twist(k,oo)
            else
              azim_chan(k1,oo)=-45.0+azim(k,oo)+twist(k,oo)
            endif
          else
            azim_chan(k1,oo)=azim(k,oo)+twist(k,oo)+shear(k,oo)
          end if
          if(azim_chan(k1,oo).gt.90.0) azim_chan(k1,oo) =
     &                              azim_chan(k1,oo)-180.0
          if(azim_chan(k1,oo).lt.-90.0) azim_chan(k1,oo) =
     &                                azim_chan(k1,oo)+180.0

c...calculate skew
          skews(jrel,k1,oo) = atand(sngl( cdabs(ztmp2(1,1)+ztmp2(2,2))/
     &                            cdabs(ztmp2(2,1)-ztmp2(1,2))) )
        
        enddo


c...write out twist and shear, and meas vs. calc impedance, for jrel=0
        if( jrel.eq.0 ) then
          if( EDI ) then
            write(*,*) ' ... ',tsite(1:index(tsite,'.edi')-1),' ...'
          else
            write(*,*) ' ... ',tsite(1:index(tsite,'.dat')-1),' ...'
          endif
          write(*,*) 'RMS Error = ', sqrt(sngl(error1)/8.)
          write(*,*) 'Shear     = ', atand(sngl(e))
          write(*,*) 'Twist     = ', atand(sngl(t))

          if( EDI ) then
            dpos = index(insite(oo),'.edi')
          else
            dpos = index(insite(oo),'.dat')
          endif
          tsite = insite(oo)

          if( FIRST_TIME ) then
            write(*,'(a)')'Opening impedance output file: '
     &                        //'imped'//tsite(4:dpos)//'dat'
            open(unit=OUTD2,file='imped'//tsite(4:dpos)//'dat',
     &         status='unknown')
            write(OUTD2,*)'     MEASURED                       '//
     &                    'ESTIMATED'
          else
            open(unit=OUTD2,file='imped'//tsite(4:dpos)//'dat',
     &           status='old',access='APPEND',err=88)
            goto 89
88          open(unit=OUTD2,file='imped'//tsite(4:dpos)//'dat',
     &         status='unknown')
89          continue
          endif

          do k = 1, count(o)
            if( SCALE ) then
              factor = sqrt(period(k,o)/(MU0*2.*PI))
            else
              factor = 1.
            endif
            write(OUTD2,*) period(k,o), factor
            write(OUTD2,2088) 
     &        real(zmeas(1,1,k,o)),imag(zmeas(1,1,k,o)),
     &        real(zmeas(1,2,k,o)),imag(zmeas(1,2,k,o)),
     &        real(zfit (1,1,k,o)),imag(zfit (1,1,k,o)),
     &        real(zfit (1,2,k,o)),imag(zfit (1,2,k,o))
            write(OUTD2,2088) 
     &        real(zmeas(2,1,k,o)),imag(zmeas(2,1,k,o)),
     &        real(zmeas(2,2,k,o)),imag(zmeas(2,2,k,o)),
     &        real(zfit (2,1,k,o)),imag(zfit (2,1,k,o)),
     &        real(zfit (2,2,k,o)),imag(zfit (2,2,k,o))
            write(OUTD2,2089) std_dev(1,1,k,o),std_dev(1,2,k,o)
            write(OUTD2,2089) std_dev(2,1,k,o),std_dev(2,2,k,o)
2088    format(4e20.8,/,4e20.8)
2089    format(2f13.7)
          enddo
          close(unit=OUTD2)

        endif

c --End of site loop

400   continue

      FIRST_TIME = .FALSE.

c --Write out strike from first run (usually main estimate)

      if( jrel.eq.0 ) then
        write(*,*) ' '
        write(*,*) 'Strike    = ', azims(0,kstart+1,1),
     &             '    (', mod(azims(0,kstart+1,1),360.),')'
        write(*,*) ' '
      endif

500   continue
c...END BIG ESTIMATE jrel LOOP 500
c...............................................................................

      kstart = kstart + bands(iband)

1000  continue
c...END BIG ESTIMATE nband LOOP 1000
c...............................................................................

      ktot = kstart
      write(*,'(a,i4)')'Total number of estimates is =', ktot

c===============================================================================
c...write main estimates (estimates 0) into azim, shear & twist arrays

      do o = 1, nsitetot
        oo = isite(o)
        do k = 1, ktot
          azim (k,oo) = azims (0,k,oo)
          shear(k,oo) = shears(0,k,oo)
          skew (k,oo) = skews (0,k,oo)
          err  (k,oo) = errs  (0,k,oo)
          twist(k,oo) = twists(0,k,oo)
          rho_a(k,oo) = rho_as(0,k,oo)
          rho_b(k,oo) = rho_bs(0,k,oo)
          pha_a(k,oo) = pha_as(0,k,oo)
          pha_b(k,oo) = pha_bs(0,k,oo)
        enddo
      enddo

      if( .not.STATS ) goto 888
      write(*,*)'Deriving errors on decomp parameters...'


c...derive standard error of azims, shears & twists, rho_a, rho_b,
c                            pha_a, pha_b

c   if no. of realizations < MINREL then take maximum & minimum values
c   *** does not include the first estimate ***

      if( nrel.lt.MINREL ) then
        write(*,*)'No. of realizations ', nrel,
     &            '  < minimum for statistics', MINREL
        write(*,*)'Extremes written as bounds to .dcmp file'
      endif

      do o = 1, nsitetot
        oo = isite(o)
        kstart = 0
        do iband = 1, nband
          do k = kstart+1, kstart+bands(iband)

c...azimuths: loop around 90 deg quadrants
          write(*,*)'...azimuth errors estimation'
          do jrel = 0, nrel
            array(jrel) = azims(jrel,k,oo) - azims(0,k,oo)
1100        if( array(jrel).ge.45. ) then
              array(jrel) = array(jrel) - 90.
              goto 1100
            else if( array(jrel).le.-45. ) then
              array(jrel) = array(jrel) + 90.
              goto 1100
            endif
            if(DEBUG) write(30,*) 'Site, freq, est, azim =', 
     &             oo,k,jrel,azims(jrel,k,o), array(jrel)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme( array, nrel-1, az_b(1,k,oo), az_b(2,k,oo) )
            az_b(2,k,oo) = azim(k,oo) + az_b(2,k,oo)
            az_b(1,k,oo) = azim(k,oo) + az_b(1,k,oo)
            if(DEBUG) write(30,*)'         azim extremes =',
     &                az_b(1,k,oo), az_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(30,*)'               azim sd =', sqrt(var)
            az_b(2,k,oo) = azim(k,oo) - sqrt(var)
            az_b(1,k,oo) = azim(k,oo) + sqrt(var)
          endif

c...shears
          write(*,*)'...shear errors estimation'
          do jrel = 0, nrel
            array(jrel) = shears(jrel,k,oo)
            if(DEBUG) write(35,*) 'Site, freq, est, shear =', 
     &                   o,k,jrel,shears(jrel,k,oo)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme( array, nrel-1, sh_b(1,k,oo), sh_b(2,k,oo) )
            if(DEBUG) write(35,*)'        shear extremes =',
     &                sh_b(1,k,oo), sh_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(35,*)'              shear sd =', sqrt(var)
            sh_b(2,k,oo) = shear(k,oo) - sqrt(var)
            sh_b(1,k,oo) = shear(k,oo) + sqrt(var)
          endif

c...twists
          write(*,*)'...twist errors estimation'
          do jrel = 0, nrel
            array(jrel) = twists(jrel,k,oo)
            if(DEBUG) write(36,*) 'Site, freq, est, twist =', 
     &                   o,k,jrel,twists(jrel,k,oo)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme( array, nrel-1, tw_b(1,k,oo), tw_b(2,k,oo) )
            if(DEBUG) write(36,*)'        twist extremes =',
     &                tw_b(1,k,oo), tw_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(36,*)'              twist sd =', sqrt(var)
            tw_b(2,k,oo) = twist(k,oo) - sqrt(var)
            tw_b(1,k,oo) = twist(k,oo) + sqrt(var)
          endif

c...errors
          write(*,*)'...error errors estimation'
          do jrel = 0, nrel
            array(jrel) = errs(jrel,k,oo)
            if(DEBUG) write(37,*) 'Site, freq, est, err =', 
     &                   o,k,jrel,errs(jrel,k,oo)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme( array, nrel-1, er_b(1,k,oo), er_b(2,k,oo) )
            if(DEBUG) write(37,*)'        error extremes =',
     &                er_b(1,k,oo), er_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(37,*)'              err sd =', sqrt(var)
            er_b(2,k,oo) = max(0.1, err(k,oo) - sqrt(var))
            er_b(1,k,oo) =          err(k,oo) + sqrt(var)
          endif

c...rho_as
          write(*,*)'...rho_a errors estimation'
          do jrel = 0, nrel
            array(jrel) = alog10(rho_as(jrel,k,oo))
            if(DEBUG) write(31,*) 'Site, freq, est, rho_a =', 
     &                   o,k,jrel,rho_as(jrel,k,oo)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme(array,nrel-1,rho_a_b(1,k,oo),rho_a_b(2,k,oo))
            rho_a_b(1,k,oo) = 10.**(rho_a_b(1,k,oo))
            rho_a_b(2,k,oo) = 10.**(rho_a_b(2,k,oo))
            if(DEBUG) write(34,*)'        rho_a extremes =',
     &                rho_a_b(1,k,oo),rho_a_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(31,*)'              rho_a sd =', sqrt(var)
            rho_a_b(2,k,oo) = 10.**(alog10(rho_a(k,oo)) - sqrt(var))
            rho_a_b(1,k,oo) = 10.**(alog10(rho_a(k,oo)) + sqrt(var))
          endif

c...rho_bs
          write(*,*)'...rho_b errors estimation'
          do jrel = 0, nrel
            array(jrel) = alog10(rho_bs(jrel,k,oo))
            if(DEBUG) write(32,*) 'Site, freq, est, rho_b =', 
     &                   o,k,jrel,rho_bs(jrel,k,oo)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme(array,nrel-1,rho_b_b(1,k,oo),rho_b_b(2,k,oo))
            rho_b_b(1,k,oo) = 10.**(rho_b_b(1,k,oo))
            rho_b_b(2,k,oo) = 10.**(rho_b_b(2,k,oo))
            if(DEBUG) write(34,*)'        rho_b extremes =',
     &                rho_b_b(1,k,oo),rho_b_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(32,*)'              rho_b sd =', sqrt(var)
            rho_b_b(2,k,oo) = 10.**(alog10(rho_b(k,oo)) - sqrt(var))
            rho_b_b(1,k,oo) = 10.**(alog10(rho_b(k,oo)) + sqrt(var))
          endif

c...pha_as
          write(*,*)'...pha_a errors estimation'
          do jrel = 0, nrel
            array(jrel) = pha_as(jrel,k,oo)
            if(DEBUG) write(33,*) 'Site, freq, est, pha_a =', 
     &                   o,k,jrel,pha_as(jrel,k,oo)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme(array,nrel-1,pha_a_b(1,k,oo),pha_a_b(2,k,oo))
            if(DEBUG) write(34,*)'        pha_a extremes =',
     &                pha_a_b(1,k,oo),pha_a_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(33,*)'              pha_a sd =', sqrt(var)
            pha_a_b(2,k,oo) = pha_a(k,oo) - sqrt(var)
            pha_a_b(1,k,oo) = pha_a(k,oo) + sqrt(var)
          endif

c...pha_bs
          write(*,*)'...pha_b errors estimation'
          do jrel = 0, nrel
            array(jrel) = pha_bs(jrel,k,oo)
            if(DEBUG) write(34,*) 'Site, freq, est, pha_b =', 
     &                   o,k,jrel,pha_bs(jrel,k,oo)
          enddo
          if( nrel.lt.MINREL ) then
            call extreme(array,nrel-1,pha_b_b(1,k,oo),pha_b_b(2,k,oo))
            if(DEBUG) write(34,*)'        pha_b extremes =',
     &                pha_b_b(1,k,oo),pha_b_b(2,k,oo)
          else
            xmean = -999.
            call jkvar( array, nrel, var )
            if(DEBUG) write(34,*)'              pha_b sd =', sqrt(var)
            pha_b_b(2,k,oo) = pha_b(k,oo) - sqrt(var)
            pha_b_b(1,k,oo) = pha_b(k,oo) + sqrt(var)
          endif

          enddo
          kstart = kstart + bands(iband)
        enddo
      enddo

c --Calculate chi-squared probabilities for a single band
      
888   if( nband.eq.1 ) then
        write(*,*)'Calculating chi^2 probability...'
        p = 0.683
        nu = neq - n
        write(*,*) ' '
        write(*,*) ' Model has ', nu, ' degrees of freedom'
        chiv = chiprob(nu,p)
        write(*,*) ' Chi-squared for 68.3% = ', chiv
        p = 0.954
        chiv = chiprob(nu,p)
        write(*,*) ' Chi-squared for 95.4% = ', chiv
        write(*,*) ' Chi-squared error = ', sngl(objf*2)
        prob = gammq(nu/2.0,sngl(objf))
        write(*,*) ' Probability of worse fit = ',100.0*prob,' %'
      endif


c===============================================================================
c...write information out to files ".dcmp"

      write(*,*)'Writing information to .dcmp files'
      do o = 1, nsitetot

c...find last directory slash (both Unix & DOS) in infile
        do i = lnblnk(infile(o)), 1, -1
          if( infile(o)(i:i).eq.fslash .or.
     &        infile(o)(i:i).eq.bslash ) then
            islash = i
            goto 58
          endif
        enddo
        islash = 0

58      tsite = infile(o)(islash+1:)
        if( EDI ) then
          dpos = index(tsite,'.edi')
        else
          dpos = index(tsite,'.dat')
        endif

c...write out parameters to ".dcmp" file

        if( DEBUG ) write(*,*)'Opening output file ',
     &                    tsite(1:dpos)//'dcmp'
        open(unit=OUTD1,file=tsite(1:dpos)//'dcmp',status='unknown')

c...open input file and write comment & info blocks to output file
c
        write(OUTD1,'(a)') '# Written by strike - version = '//ver
        write(OUTD1,'(a)') '# date: '//date_string
        write(OUTD1,'(a)') '#'
        write(OUTD1,'(a)') '# input file >'//tsite(1:lnblnk(tsite))
        write(OUTD1,'(a)') '#'
        write(OUTD1,'(a)') '# Normalization type : '//normtype
        write(OUTD1,'(a)') '#'
        write(OUTD1,'(a)') '# Bands and objective functions achieved:'
        write(OUTD1,'(a)') '# Band     Per-min      Per-max'//
     &                 '             Start           End'
        do iband = 1, nband
          write(OUTD1,'(a,i4,2g15.4,2f15.1)') '#', iband,  
     &                  perband(1,iband), perband(2,iband),
     &                  objstart(iband), objend(iband)
        enddo
        write(OUTD1,'(a)') '#'
        write(OUTD1,'(a)') '# bounds on distortion parameters:'
        write(OUTD1,'(2(a,f8.1))') '# azimuth:',
     $       sngl(alow), ' - ',sngl(ahigh)
        write(OUTD1,'(2(a,f8.1))') '# shear  :',
     $       sngl(slow), ' - ', sngl(shigh)
        write(OUTD1,'(2(a,f8.1))') '# twist  :',
     $       sngl(tlow), ' - ', sngl(thigh)
        write(OUTD1,'(a)') '#'
        call trunc(inlist,linsite)
        write(OUTD1,'(a,a)') '# Site list : ',inlist(1:linsite)
        if( .not.SINGLE ) then
          do j2 = 1, nsite, 5
            write(OUTD1,'(a,5(a12,1x))') '#             ', (insite(j),
     $          j=j2,min0(j2+4,nsite))
          enddo
        endif
        write(OUTD1,'(a)') '#'


c...write comment and information blocks to output file
        if( EDI ) then
          close(unit=IND1)
          open (unit=IND1,file=infile(o),status='old',err=59,iostat=ios)
          call edicopycomm(IND1,OUTD1)
          call edicopyinfo(IND1,OUTD1)
          close(unit=IND1)
          goto 60
        else
c...reopen input file to copy comment and info blocks to output file...
          close(unit=IND1)
          open (unit=IND1,file=infile(o),status='old',err=59,iostat=ios)
          call jcopycomm(IND1,OUTD1)
          call jcopyinfo(IND1,OUTD1)
          close(unit=IND1)
          goto 60
        endif

59      write(*,*) 'Error opening ',
     &          filein(o)(1:lnblnk(filein(o))), '   iostat=',ios
        write(*,*)'No comment or information blocks written to '//
     &            'output file'



c...derive no. of estimates for this site
60      jtot = 0
        j1 = 0
        j2 = 0
        do j = 1, ktot
          if( azim(j,o).ne.-9999. ) then
            jtot = jtot+1
            if( j1.eq.0 ) j1 = j
            j2 = j
          endif
        enddo

        write(*,'(a,i3)')'Writing dcmp info to file...'//
     &           tsite(1:dpos)//'dcmp    no. estimates=',jtot

c...azimuth
        write(OUTD1,*)3
        write(OUTD1,*) jtot,' regional azimuth (cf AZIMUTH coordinate',
     &                   ' frame)'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),mod(azim(j,o),360.),
     &                mod(az_b(1,j,o),360.),
     &                mod(az_b(2,j,o),360.)
        enddo

c...shear
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' shear angle'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),shear(j,o),
     &                sh_b(1,j,o),sh_b(2,j,o)
        enddo

c...channelling angle
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' channelling angle'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o), azim_chan(j,o),
     &                chan_b(1,j,o),chan_b(2,j,o)
        enddo

c...twist
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' twist angle'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),twist(j,o),
     &                tw_b(1,j,o),tw_b(2,j,o)
        enddo

c...apparent resistivity a
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' app rho a'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),rho_a(j,o),
     &                rho_a_b(1,j,o),rho_a_b(2,j,o)
        enddo

c...apparent resistivity b
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' app rho b'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),rho_b(j,o),
     &                rho_b_b(1,j,o),rho_b_b(2,j,o)
        enddo

c...impedance phase a
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' imped phase a'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),pha_a(j,o),
     &                pha_a_b(1,j,o),pha_a_b(2,j,o)
        enddo

c...impedance phase b
        write(OUTD1,*)3
        write(OUTD1,*)ktot,' imped phase b'
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),pha_b(j,o),
     &                pha_b_b(1,j,o),pha_b_b(2,j,o)
        enddo

c...error normalize error to rms and (calculate total rms error)
        toterr = 0.
        do j = j1, j2
           err(j,o) = sqrt(err(j,o)/8.)
           er_b(1,j,o) = sqrt(er_b(1,j,o)/8.)
           er_b(2,j,o) = sqrt(er_b(2,j,o)/8.)
           toterr = toterr + err(j,o)
        enddo
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' av. rms error ', toterr/float(j2-j1+1)
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),err(j,o),
     &                er_b(1,j,o),er_b(2,j,o)
        enddo    

c...skew
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' skew '
        do j = j1, j2
           write(OUTD1,*)freqs(j,o),skew(j,o),
     &                sk_b(1,j,o),sk_b(2,j,o)
        enddo

c...anisotropy
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' anis '
        do j = j1, j2
           if( rho_b(j,o).ne.0. ) then
             zrat = sqrt( rho_a(j,o)/rho_b(j,o) )
             s  = (zrat-1.)/(zrat+1.)
             if( rho_b_b(2,j,o).ne.0. ) then
               zrat = sqrt(rho_a_b(1,j,o)/rho_b_b(2,j,o) )
               sl = (zrat-1.)/(zrat+1.)
             else
               sl = 0.
             endif
             if( rho_b_b(1,j,o).ne.0. ) then
               zrat = sqrt(rho_a_b(2,j,o)/rho_b_b(1,j,o) )
               su = (zrat-1.)/(zrat+1.)
             else
               su = 0.
             endif
           else
             s  = 0.
             sl = 0.
             su = 0.
           endif
           write(OUTD1,*)freqs(j,o),s, sl, su
        enddo

c...phase difference
        write(OUTD1,*)3
        write(OUTD1,*)jtot,' phadif '
        do j = j1, j2
           s = pha_a(j,o) - pha_b(j,o)
           sl = pha_a_b(1,j,o) - pha_b_b(2,j,o)
           su = pha_a_b(2,j,o) - pha_b_b(1,j,o)
           write(OUTD1,*)freqs(j,o),s, sl, su
        enddo

      enddo

      stop

9998  write(*,*)'Error opening ',inlist(1:lnblnk(inlist))
      stop

9999  write(*,*)'Error opening ',infile(o)(1:lnblnk(infile(o)))
      stop
      end
