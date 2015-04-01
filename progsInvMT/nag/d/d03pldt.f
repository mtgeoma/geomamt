      SUBROUTINE D03PLD(NPDE,T,X,ULEFT,URIGHT,RFLUX,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Dummy routine used as an external in D03PFR, used only in D03PLF
C  and D03PSF.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  RFLUX(NPDE), ULEFT(NPDE), URIGHT(NPDE)
C     .. Executable Statements ..
      RETURN
      END
