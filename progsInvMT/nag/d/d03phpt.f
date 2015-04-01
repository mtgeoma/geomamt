      SUBROUTINE D03PHP(NPDE,T,X,U,DUDX,C,F,R,IRES)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C----------------------------------------------------------------------
C  A dummy routine as an external in D03PCH used only for D03PHF
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE,NPDE), DUDX(NPDE), F(NPDE), R(NPDE),
     *                  U(NPDE)
C     .. Executable Statements ..
      RETURN
      END
