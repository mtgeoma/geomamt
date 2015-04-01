      SUBROUTINE D03PKN(NPDE,T,IBND,NOBC,U,UDOT,NV,V,VDOT,RES,IRES)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C    -------------------------------------------------------------------
C     A dummy routine as an external in D03PEG used only for D03PEF
C    -------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IBND, IRES, NOBC, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  RES(NPDE), U(NPDE), UDOT(NPDE), V(*), VDOT(*)
C     .. Executable Statements ..
      RETURN
      END
