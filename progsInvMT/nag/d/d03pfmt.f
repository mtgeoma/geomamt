      SUBROUTINE D03PFM(NPDE,T,X,U,UX,P,C,D,S,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C----------------------------------------------------------------------
C  A dummy routine used as an external in D03PFH used only for
C  D03PLF and D03PSF
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE), D(NPDE), P(NPDE,NPDE), S(NPDE),
     *                  U(NPDE), UX(NPDE)
C     .. Executable Statements ..
      RETURN
      END
