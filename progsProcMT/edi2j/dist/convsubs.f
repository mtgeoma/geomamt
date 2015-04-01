c-----------------------------------------------------------------------
      complex function convr2z( rho, pha, per )

      real PI, MU
      parameter( PI = 3.141592 , MU = 4.*PI*1.e-7 )

      real cosd, sind

c---converts rho & pha to Z in S.I. 
c      units of Z are ohms (E/H)

      omega = 2.*PI/per

      zm = sqrt(omega*rho*MU)
      convr2z  = cmplx(zm*cosd(pha),zm*sind(pha))
      return
      end

c-----------------------------------------------------------------------
      real function convz2r( z, per )
      complex z

      real PI, MU
      parameter( PI = 3.141592 , MU = 4.*PI*1.e-7 )

      omega = 2.*PI/per

      convz2r  = cabs(z)**2. / (omega * MU)
      return
      end

c-----------------------------------------------------------------------
      real function convz2p( z, per )
      complex z

      real PI, MU
      parameter( PI = 3.141592 , MU = 4.*PI*1.e-7 )

      real atan2d

      omega = 2.*PI/per

      convz2p = atan2d( aimag(z), real(z) )
      return
      end

c-----------------------------------------------------------------------
      complex function convr2c( rho, pha, per )

      real PI, MU
      complex IIMAG
      parameter ( 
     &            PI   = 3.141592 , 
     &            MU   = 4.*PI*1.e-7,
     &            IIMAG = ( 0., 1.)
     &           )

      omega = 2.*PI/per

      zm = sqrt(omega*rho/MU)
      convr2c  = cmplx(zm*cosd(pha),zm*sind(pha))/(IIMAG*omega)
      return
      end
