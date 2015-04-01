      SUBROUTINE D03PEH(PKFPDE,PDEPEF,T,X,NPDE,U,UDOT,DUDX,NV,V,VDOT,
     *                  RES,IRES,IIFLAG)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     ----------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     PDE description routine to call external routines in D03PKZ
C     used in D03PEF, D03PKF.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IIFLAG, IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  DUDX(NPDE), RES(NPDE), U(NPDE), UDOT(NPDE),
     *                  V(*), VDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          PDEPEF, PKFPDE
C     .. Executable Statements ..
      CALL PKFPDE(NPDE,T,X,U,UDOT,DUDX,RES,IRES)
      CALL PDEPEF(NPDE,T,X,U,UDOT,DUDX,NV,V,VDOT,RES,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
      END
