      SUBROUTINE D03PCG(PHFBND,BNDPCF,T,BETA,GAMMA,U,UX,NPDE,I,NV,V,VD,
     *                  IRES,IIFLAG)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  This subroutine calls the external routines provided by the user
C  of D03PCF or D03PHF to describe the boundary conditions.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           I, IIFLAG, IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  BETA(NPDE), GAMMA(NPDE), U(NPDE), UX(NPDE),
     *                  V(*), VD(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPCF, PHFBND
C     .. Executable Statements ..
      CALL PHFBND(NPDE,T,U,UX,I,BETA,GAMMA,IRES)
      CALL BNDPCF(NPDE,T,U,UX,NV,V,VD,I,BETA,GAMMA,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
      END
