      SUBROUTINE D03PJJ(NPDE,T,U,UX,IB,BETA,GAMMA,IRES)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ---------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine calls the fixed name routine provided by the user
C     of D03PDF to describe the boundary conditions.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ---------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IB, IRES, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  BETA(NPDE), GAMMA(NPDE), U(NPDE), UX(NPDE)
C     .. Executable Statements ..
      RETURN
C
      END
