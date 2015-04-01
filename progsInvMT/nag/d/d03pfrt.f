      SUBROUTINE D03PFR(FLXPFF,FLXPLF,NPDE,T,X,NV,V,ULEFT,URIGHT,RFLUX,
     *                  IRES,IIFLAG)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Routine to call user-supplied flux routine.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IIFLAG, IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  RFLUX(NPDE), ULEFT(NPDE), URIGHT(NPDE), V(*)
C     .. Subroutine Arguments ..
      EXTERNAL          FLXPFF, FLXPLF
C     .. Executable Statements ..
      CALL FLXPFF(NPDE,T,X,ULEFT,URIGHT,RFLUX,IRES)
      CALL FLXPLF(NPDE,T,X,NV,V,ULEFT,URIGHT,RFLUX,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
      END
