      SUBROUTINE D03PKT(PKFPDE,PDEPEF,T,X,NPDE,U,DUDX,C,F,R,NV,V,VDOT,
     *                  IRES,IIFLAG)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IIFLAG, IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE,NPDE), DUDX(NPDE), F(NPDE), R(NPDE),
     *                  U(NPDE), V(*), VDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          PDEPEF, PKFPDE
C     .. Local Scalars ..
      INTEGER           I
C     .. Executable Statements ..
      DO 20 I = 1, NPDE
         R(I) = 0.0D0
   20 CONTINUE
      RETURN
      END
