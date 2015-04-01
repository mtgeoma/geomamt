      subroutine edizdatin( ind, n, per, zxx, zxxe, zxy, zxye,
     &                                   zyx, zyxe, zyy, zyye )

c---------------------------------------------------------------------

      implicit none

c---------------------------------------------------------------------

      include 'ctrlblk.inc'
      include 'size.inc'

c---------------------------------------------------------------------

      real PI
      parameter ( PI = 3.141592 )
      real MU
      parameter ( MU = 4.*PI*1.E-7 )

c---------------------------------------------------------------------

      integer ind, n
      complex zxx(*), zxy(*), zyx(*), zyy(*)
      real per(*), zxxe(*), zxye(*), zyxe(*), zyye(*)

c---------------------------------------------------------------------

      COMPLEX Z(3,2), zdata(6,MAXDAT)
      real zvar(6,MAXDAT)
      REAL rho(2,2), phs(2,2), srh(2,2), sph(2,2)
      REAL vimp(3,2), delta, avgs, angle, freq(MAXDAT)
      REAL sdata(7,7,MAXDAT), info(MAXDAT,3), temp(7,7)
      real sdzxx, sdzxy, sdzyx, sdzyy
      integer i, j, k

      character*80 line
      character*6 site

      logical SPECTRA
      logical DEBUG
      logical INC
	
      real factor

      common /cmndbg/ DEBUG

c-------------------------------------------------------------------------------
c...initialize

c...factor to convert impedances from field units to S.I.
      factor = 4.*PI*1.e-4

      angle = 0.0

c-------------------------------------------------------------------------------
c...check whether MTPARAMS file or SPECTRA file by looking for >FREQ section

      if( DEBUG ) write(*,*)'edizdatin: entered'

      read(ind,'(a)',end=10) line
      do while(index(line,'>FREQ ').eq.0)
        read(ind,'(a)',end=10) line
      enddo
      if( index(line,'ORDER=INC').ne.0 ) then
        INC = .TRUE.
      else
        INC = .FALSE.
      endif


c...read in freqs
      if( DEBUG ) write(*,*)'edizdatin: MTPARAMS edi file to be read'
      SPECTRA = .FALSE.
      goto 20

c...SPECTRA file. Have to get freqs off each SPECTRA section
10    if( DEBUG ) write(*,*)'edizdatin: SPECTRA edi file to be read'
      SPECTRA = .TRUE.


c-------------------------------------------------------------------------------
c...read edi-file

20    rewind( ind )
      if( SPECTRA ) then
        call edispecread(ind,freq,n,sdata,info,site)
        if( DEBUG ) then
          write(*,*)'edizdatin: return from edispecread routine: n =', n
        endif
      else
        call ediresread(ind,freq,n,zdata,zvar,info,site)
        if( DEBUG ) then
          write(*,*)'edizdatin: return from ediresread routine: n =', n
        endif
      endif

c-------------------------------------------------------------------------------
c-- loop over frequencies 

      do i = 1, n

        if( SPECTRA ) then
          do j = 1, 7
            do k = 1, 7
              temp(j,k) = sdata(j,k,i)
            enddo
          enddo
          avgs =  info(i,3)
          delta = angle - info(i,1)
          call mtcomp(temp,freq(i),delta,avgs,Z,vimp,rho,phs,
     &                srh,sph)
          sdzxx = sqrt(vimp(1,1))
          sdzxy = sqrt(vimp(1,2))
          sdzyx = sqrt(vimp(2,1))
          sdzyy = sqrt(vimp(2,2))

        else
c...convert impedances from field units to S.I.
          if( factor.ne.1. ) then
            Z(1,1) = factor * zdata(1,i)
            Z(1,2) = factor * zdata(2,i)
            Z(2,1) = factor * zdata(3,i)
            Z(2,2) = factor * zdata(4,i)
            sdzxx = factor * sqrt(zvar(1,i))
            sdzxy = factor * sqrt(zvar(2,i))
            sdzyx = factor * sqrt(zvar(3,i))
            sdzyy = factor * sqrt(zvar(4,i))
          else
            sdzxx = sqrt(zvar(1,i))
            sdzxy = sqrt(zvar(2,i))
            sdzyx = sqrt(zvar(3,i))
            sdzyy = sqrt(zvar(4,i))
          endif

        endif

        zxx(i) = Z(1,1)
        zxy(i) = Z(1,2)
        zyx(i) = Z(2,1)
        zyy(i) = Z(2,2)

        zxxe(i) = sdzxx
        zxye(i) = sdzxy
        zyxe(i) = sdzyx
        zyye(i) = sdzyy

        per(i) = 1./freq(i)

      enddo

c...reverse order if increasing frequency
      if( INC ) then
        call  revers( n, per )
        call crevers( n, zxx )
        call crevers( n, zxy )
        call crevers( n, zyx )
        call crevers( n, zyy )
        call  revers( n, zxxe )
        call  revers( n, zxye )
        call  revers( n, zyxe )
        call  revers( n, zyye )
      endif


      if( DEBUG ) write(*,*)'edizdatin: exiting'

      return
      end




           

