      SUBROUTINE D03PEG(PKFBND,BNDPEF,T,IBND,NPDE,U,UDOT,NV,V,VDOT,NOBC,
     *                  RES,IRES,IIFLAG)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     ----------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine calls the external routines provided by the
C     user of D03PEF or D03PKF to describe the boundary conditions.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IBND, IIFLAG, IRES, NOBC, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  RES(NOBC), U(NPDE), UDOT(NPDE), V(*), VDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPEF, PKFBND
C     .. Executable Statements ..
      CALL PKFBND(NPDE,T,IBND,NOBC,U,UDOT,RES,IRES)
      CALL BNDPEF(NPDE,T,IBND,NOBC,U,UDOT,NV,V,VDOT,RES,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
      END
