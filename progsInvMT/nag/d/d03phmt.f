      SUBROUTINE D03PHM(NPDE,T,X,U,DUDX,NV,V,VDOT,C,F,R,IRES)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C----------------------------------------------------------------------
C  A dummy routine as an external in D03PCH used only for D03PCF
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE,NPDE), DUDX(NPDE), F(NPDE), R(NPDE),
     *                  U(NPDE), V(*), VDOT(*)
C     .. Executable Statements ..
      RETURN
      END
