      SUBROUTINE D03PKM(NPDE,T,X,U,UDOT,DUDX,NV,V,VDOT,RES,IRES)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C ----------------------------------------------------------------------
C     A dummy routine as an external in D03PEH used only for D03PEF
C ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  DUDX(NPDE), RES(NPDE), U(NPDE), UDOT(NPDE),
     *                  V(*), VDOT(*)
C     .. Executable Statements ..
      RETURN
      END
