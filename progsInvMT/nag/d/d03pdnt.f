      SUBROUTINE D03PDN(PDEFCN,BNDARY,DUMPD1,DUMPD2,NEQN,T,HLAST,H,Y,
     *                  YDOT,YSAVE,NYH,R,ACOR,RESWK,NRESWK,WKMON,NWKMON,
     *                  IMON,INLN,HMIN,HMAX,EVNRMF,PDEFN,BNDR,IRES,
     *                  IXFIX,NIXFIX)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C ----------------------------------------------------------------------
C     Dummy monitor routine for Chebyshev discretisation.
C ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HLAST, HMAX, HMIN, T
      INTEGER           IMON, INLN, IRES, NEQN, NIXFIX, NRESWK, NWKMON,
     *                  NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(*), R(*), RESWK(NRESWK), WKMON(NWKMON),
     *                  Y(*), YDOT(*), YSAVE(NYH,*)
      INTEGER           IXFIX(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, BNDR, DUMPD1, DUMPD2, EVNRMF, PDEFCN,
     *                  PDEFN
C     .. Executable Statements ..
      RETURN
C
      END
