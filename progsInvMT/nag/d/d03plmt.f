      SUBROUTINE D03PLM(NPDE,T,X,U,UX,NV,V,VDOT,P,C,D,S,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C----------------------------------------------------------------------
C  A dummy routine used as an external in D03PFH used only for
C  D03PFF.
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE), D(NPDE), P(NPDE,NPDE), S(NPDE),
     *                  U(NPDE), UX(NPDE), V(*), VDOT(*)
C     .. Executable Statements ..
      RETURN
      END
