      SUBROUTINE D03PDH(BND11,BND22,T,BETA,GAMMA,U,UX,NPDE,LEFT,NV,V,
     *                  VDOT,IRES,IIFLAG)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine calls the fixed name routine provided by the user
C     of D03PDF to describe the boundary conditions
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IIFLAG, IRES, NPDE, NV
      LOGICAL           LEFT
C     .. Array Arguments ..
      DOUBLE PRECISION  BETA(NPDE), GAMMA(NPDE), U(NPDE), UX(NPDE),
     *                  V(*), VDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BND11, BND22
C     .. Local Scalars ..
      INTEGER           IB
C     .. Executable Statements ..
      IF (LEFT) THEN
         IB = 0
      ELSE
         IB = 11
      END IF
C
      CALL BND11(NPDE,T,U,UX,IB,BETA,GAMMA,IRES)
      CALL BND22(NPDE,T,U,UX,NV,V,VDOT,IB,BETA,GAMMA,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
C
      END
