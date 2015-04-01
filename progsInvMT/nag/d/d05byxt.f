      SUBROUTINE D05BYX(GAMFUN,IS)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ------------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine uses the NAG function  S14AAF  to generates the
C     quotient containing the Gamma function. The values are stored
C     in  GAMFUN.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IS
C     .. Array Arguments ..
      DOUBLE PRECISION  GAMFUN(IS)
C     .. Local Scalars ..
      DOUBLE PRECISION  ARGMNT, GAMDEN, GAMNUM
      INTEGER           IFAIL, M
C     .. External Functions ..
      DOUBLE PRECISION  S14AAF
      EXTERNAL          S14AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IFAIL = 0
C
      DO 20 M = 0, IS - 1
         ARGMNT = 1.0D0 + 0.5D0*DBLE(M)
         GAMNUM = S14AAF(ARGMNT,IFAIL)
         GAMDEN = S14AAF(ARGMNT+0.5D0,IFAIL)
         GAMFUN(M+1) = GAMNUM/GAMDEN
   20 CONTINUE
C
      RETURN
      END
