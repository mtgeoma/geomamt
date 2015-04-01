      SUBROUTINE D03PFH(PLFPDE,PDEPFF,T,X,NPDE,U,UX,P,C,D,S,NV,V,VDOT,
     *                  IRES,IIFLAG)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PDE description routine to call external routines used in D03PFF,
C  D03PLF, D03PSF.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IIFLAG, IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE), D(NPDE), P(NPDE,NPDE), S(NPDE),
     *                  U(NPDE), UX(NPDE), V(*), VDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          PDEPFF, PLFPDE
C     .. Executable Statements ..
      CALL PLFPDE(NPDE,T,X,U,UX,P,C,D,S,IRES)
      CALL PDEPFF(NPDE,T,X,U,UX,NV,V,VDOT,P,C,D,S,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
      END
