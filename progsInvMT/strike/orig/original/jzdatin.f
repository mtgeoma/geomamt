      subroutine jzdatin( ind, n, per, zxx, zxxe, zxy, zxye, 
     &                                 zyx, zyxe, zyy, zyye )

c---------------------------------------------------------------------
c...modification history
c
c   agj: 7/Feb/2001:   set negative errors to zero
c   agj: 21/July/2004: checks for S.I. or FIELD or neither

c---------------------------------------------------------------------

      implicit none

c---------------------------------------------------------------------

      real PI, MU
      parameter( PI = 3.141592 , MU = 4.*PI*1.e-7 )

c---------------------------------------------------------------------

      integer ind, n
      complex zxx(*), zxy(*), zyx(*), zyy(*)
      real per(*), zxxe(*), zxye(*), zyxe(*), zyye(*)

c---------------------------------------------------------------------
      real convf2si

      character line*80
      real zr, zi, zm, rho, pha, rhomx, rhomn, phamx, phamn,
     &     err_r, err_p, zmx, zmn, period, omega
      integer i
      logical DEBUG, SI

      real cosd, sind, tand

      common /cmndbg/ DEBUG

c---------------------------------------------------------------------
c...initialize

      convf2si = sqrt(0.2 * (2*PI) * (4*PI*1.e-7) )

c---------------------------------------------------------------------

      if( DEBUG ) write(*,*)'jzdatin: entered'

c...reads in a J-format file of impedances
c   if impedances (ZXX, ZXY etc.) not found, then reads in apparent
c   resistivities and phases and converts them

c...search for ZXX in file
      read(ind,'(a)') line
      do while(line(1:3).ne.'ZXX')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'ZXX found - impedances to be read in'

c...check whether units designation present (S.I. or SI) or (FIELD or field)
      if( index(line,'S.I.').ne.0 .or. 
     &    index(line,'SI').ne.0 ) then
        SI = .TRUE.
      elseif( index(line,'FIELD').ne.0 .or. 
     &        index(line,'field').ne.0 ) then
        SI = .FALSE.
      else
        write(*,*)'Ambiguous units for impedance'
        SI = .TRUE.
        call lin('S.I. units', SI )
      endif
      if(DEBUG)write(*,'(a)')'Impedance units in SI =', SI


      read(4,*) n
      do i = 1, n
        read(4,*) period, zr, zi, zxxe(i)
        zxx(i) = cmplx(zr, zi)
        if( zxxe(i).lt.0. ) zxxe(i) = 0.
      enddo
      rewind(unit=ind)

      read(ind,'(a)') line
      do while(line(1:3).ne.'ZXY')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'ZXY found - impedances to be read in'
      read(4,*) n
      do i = 1, n
        read(4,*) period, zr, zi, zxye(i)
        zxy(i) = cmplx(zr, zi)
        if( zxye(i).lt.0. ) zxye(i) = 0.
      enddo
      rewind(unit=ind)

      read(ind,'(a)') line
      do while(line(1:3).ne.'ZYX')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'ZYX found - impedances to be read in'
      read(4,*) n
      do i = 1, n
        read(4,*) period, zr, zi, zyxe(i)
        zyx(i) = cmplx(zr, zi)
        if( zyxe(i).lt.0. ) zyxe(i) = 0.
      enddo
      rewind(unit=ind)

      read(ind,'(a)') line
      do while(line(1:3).ne.'ZYY')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'ZYX found - impedances to be read in'
      read(4,*) n
      do i = 1, n
        read(4,*) period, zr, zi, zyye(i)
        if( period.gt.0. ) then
          per(i) = period
        else
          per(i) = -1./period
        endif
        zyy(i) = cmplx(zr, zi)
        if( zyye(i).lt.0. ) zyye(i) = 0.
      enddo

      if( .not.SI ) then
        do i = 1, n
          zxx(i)  = convf2si*zxx(i)
          zxxe(i) = convf2si*zxxe(i)
          zxy(i)  = convf2si*zxy(i)
          zxye(i) = convf2si*zxye(i)
          zyx(i)  = convf2si*zyx(i)
          zyxe(i) = convf2si*zyxe(i)
          zyy(i)  = convf2si*zyy(i)
          zyye(i) = convf2si*zyye(i)
        enddo
      endif

      return

100   rewind(unit=ind)
      if(DEBUG)write(*,'(a)')'ZXX not found - reading in rhos and phas'
      
      read(ind,'(a)') line
      do while(line(1:3).ne.'RXX')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'RXX found - rhos/phas to be read in'
      read(4,*) n
      do i = 1, n
        read(4,*) period, rho, pha, rhomx, rhomn, phamx, phamn
        if( period.gt.0. ) then
          omega = 2.*PI/period
        else
          omega = -2.*PI*period
        endif
        zm = sqrt(omega*rho*MU)
        zxx(i) = cmplx(zm*cosd(pha),zm*sind(pha))
        zmx = sqrt(omega*rhomx*MU)
        zmn = sqrt(omega*rhomn*MU)
        err_r = (zmx-zmn)/2.
        err_p = zm*tand( (phamx-phamn)/2. )
        zxxe(i) = (err_r + err_p)/2.
      enddo
      rewind(unit=ind)

      read(ind,'(a)') line
      do while(line(1:3).ne.'RXY')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'RXY found - rhos/phas to be read in'
      read(4,*) n
      do i = 1, n
        read(4,*) period, rho, pha, rhomx, rhomn, phamx, phamn
        if( period.gt.0. ) then
          omega = 2.*PI/period
        else
          omega = -2.*PI*period
        endif
        zm = sqrt(omega*rho*MU)
        zxy(i) = cmplx(zm*cosd(pha),zm*sind(pha))
        zmx = sqrt(omega*rhomx*MU)
        zmn = sqrt(omega*rhomn*MU)
        err_r = (zmx-zmn)/2.
        err_p = zm*tand( (phamx-phamn)/2. )
        zxye(i) = (err_r + err_p)/2.
      enddo
      rewind(unit=ind)

      read(ind,'(a)') line
      do while(line(1:3).ne.'RYX')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'RYX found - rhos/phas to be read in'
      read(4,*) n
      do i = 1, n
        read(4,*) period, rho, pha, rhomx, rhomn, phamx, phamn
        if( period.gt.0. ) then
          omega = 2.*PI/period
        else
          omega = -2.*PI*period
        endif
        zm = sqrt(omega*rho*MU)
        zyx(i) = cmplx(zm*cosd(pha),zm*sind(pha))
        zmx = sqrt(omega*rhomx*MU)
        zmn = sqrt(omega*rhomn*MU)
        err_r = (zmx-zmn)/2.
        err_p = zm*tand( (phamx-phamn)/2. )
        zyxe(i) = (err_r + err_p)/2.
      enddo
      rewind(unit=ind)

      read(ind,'(a)') line
      do while(line(1:3).ne.'RYY')
        read(ind,'(a)',end=100) line
      enddo
      if(DEBUG)write(*,'(a)')'RYY found - rhos/phas to be read in'
      read(4,*) n
      do i = 1, n
        read(4,*) period, rho, pha, rhomx, rhomn, phamx, phamn
        if( period.gt.0. ) then
          per(i) = period
        else
          per(i) = -1./period
        endif
        if( period.gt.0. ) then
          omega = 2.*PI/period
        else
          omega = -2.*PI*period
        endif
        zm = sqrt(omega*rho*MU)
        zyy(i) = cmplx(zm*cosd(pha),zm*sind(pha))
        zmx = sqrt(omega*rhomx*MU)
        zmn = sqrt(omega*rhomn*MU)
        err_r = (zmx-zmn)/2.
        err_p = zm*tand( (phamx-phamn)/2. )
        zyye(i) = (err_r + err_p)/2.
      enddo

      if( DEBUG ) write(*,*)'jzdatin: completed'

      return

      end
