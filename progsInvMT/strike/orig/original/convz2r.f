c*****************************************************************
      real function convz2r( z, per )
c
c   convert Z to apparent resistivity
c
c*****************************************************************
      implicit none
      complex z
      real factor,mu0,omega,per,pi
      
      pi = 3.14159265
      mu0 = 4.0*PI/(10.**7)
      omega = 2.*PI/per

c...this factor has to be used for Z in field units
c      factor=(4.0*pi)/(10.**4)
      factor = 1.0
 
      convz2r  = cabs(factor*z)**2./(mu0*omega)

      return
      end
