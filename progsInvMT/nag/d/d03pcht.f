      SUBROUTINE D03PCH(PHFPDE,PDEPCF,T,X,NPDE,U,DUDX,C,F,R,NV,V,VDOT,
     *                  IRES,IIFLAG)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PDE description routine to call external routines in D03PHZ used
C  in D03PCF, D03PHF.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IIFLAG, IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE,NPDE), DUDX(NPDE), F(NPDE), R(NPDE),
     *                  U(NPDE), V(*), VDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          PDEPCF, PHFPDE
C     .. Executable Statements ..
      CALL PHFPDE(NPDE,T,X,U,DUDX,C,F,R,IRES)
      CALL PDEPCF(NPDE,T,X,U,DUDX,NV,V,VDOT,C,F,R,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
      END
