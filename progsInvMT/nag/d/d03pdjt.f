      SUBROUTINE D03PDJ(PDE11,PDE22,T,X,NPTL,NPDE,U,DUDX,C,F,R,NV,V,
     *                  VDOT,IRES,IIFLAG)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     PDE description routine to call fixed name routine defined
C     by the user of D03PDF.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IIFLAG, IRES, NPDE, NPTL, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE,NPDE,NPTL), DUDX(NPDE,NPTL),
     *                  F(NPDE,NPTL), R(NPDE,NPTL), U(NPDE,NPTL), V(*),
     *                  VDOT(*), X(NPTL)
C     .. Subroutine Arguments ..
      EXTERNAL          PDE11, PDE22
C     .. Executable Statements ..
      CALL PDE11(NPDE,T,X,NPTL,U,DUDX,C,F,R,IRES)
      CALL PDE22(NPDE,T,X,NPTL,U,DUDX,NV,V,VDOT,C,F,R,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
C
      END
