      SUBROUTINE D03PLN(NPDE,NPTS,T,X,U,NV,V,VD,I,G,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C----------------------------------------------------------------------
C  A dummy routine used as an external in D03PFG used only for
C  D03PFF.
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           I, IRES, NPDE, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  G(NPDE), U(NPDE,NPTS), V(*), VD(*), X(NPTS)
C     .. Executable Statements ..
      RETURN
      END
