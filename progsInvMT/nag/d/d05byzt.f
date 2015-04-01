      SUBROUTINE D05BYZ(WT,IS,LENP,LIQ,SW,LSW,AUXARY,A)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ----------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++
C     D05BYZ  calculates the values of the starting
C     fractional weights associated with a BDF fractional
C     rule.
C     On exit, SW  contains the values of the
C     starting weights.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ----------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IS, LENP, LIQ, LSW
C     .. Array Arguments ..
      DOUBLE PRECISION  A(2*LIQ), AUXARY(2*LIQ), SW(LSW,0:IS-1), WT(LIQ)
C     .. Local Arrays ..
      DOUBLE PRECISION  GAMFUN(11), TEMPST(11), VANMAT(11)
C     .. External Subroutines ..
      EXTERNAL          D05BYK, D05BYY
C     .. Executable Statements ..
C
      CALL D05BYY(WT,IS,LENP,LIQ,SW,LSW,GAMFUN,AUXARY,A)
C
      CALL D05BYK(LENP,IS,SW,LSW,GAMFUN,TEMPST,VANMAT)
C
      RETURN
      END
