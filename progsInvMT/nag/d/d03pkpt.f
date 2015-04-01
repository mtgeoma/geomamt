      SUBROUTINE D03PKP(NPDE,T,X,U,UDOT,DUDX,RES,IRES)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C ----------------------------------------------------------------------
C     A dummy routine as an external in D03PEH used only for D03PKF
C ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  DUDX(NPDE), RES(NPDE), U(NPDE), UDOT(NPDE)
C     .. Executable Statements ..
      RETURN
      END
