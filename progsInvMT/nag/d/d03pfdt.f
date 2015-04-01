      SUBROUTINE D03PFD(NPDE,T,X,NV,V,ULEFT,URIGHT,RFLUX,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Dummy routine used as an external in D03PFR, used only in D03PFF.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  RFLUX(NPDE), ULEFT(NPDE), URIGHT(NPDE), V(*)
C     .. Executable Statements ..
      RETURN
      END
