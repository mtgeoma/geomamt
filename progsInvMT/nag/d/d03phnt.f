      SUBROUTINE D03PHN(NPDE,T,U,UX,NV,V,VD,I,BETA,GAMMA,IRES)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-------------------------------------------------------------------
C  A dummy routine as an external in D03PCG used only for D03PCF
C-------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           I, IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  BETA(NPDE), GAMMA(NPDE), U(NPDE), UX(NPDE),
     *                  V(*), VD(*)
C     .. Executable Statements ..
      RETURN
      END
