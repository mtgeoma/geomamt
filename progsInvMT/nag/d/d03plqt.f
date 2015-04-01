      SUBROUTINE D03PLQ(NPDE,NPTS,T,X,U,I,G,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C----------------------------------------------------------------------
C  A dummy routine used as an external in D03PFG used only for
C  D03PLF and D03PSF
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           I, IRES, NPDE, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  G(NPDE), U(NPDE,3), X(NPTS)
C     .. Executable Statements ..
      RETURN
      END
