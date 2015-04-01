      SUBROUTINE G02HKZ(EPS,NVAR,A2,B2)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES THE CONSTANTS A2 AND B2 FOR THE
C     HUBER WEIGHT FUNCTION USED IN ROBUST
C     COVARIANCE ESTIMATION
C
C     .. Parameters ..
      DOUBLE PRECISION  TOL
      PARAMETER         (TOL=0.00005D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A2, B2, EPS
      INTEGER           NVAR
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, EP, FA, FB, TOLX, XP
      INTEGER           IFAULT, IND
C     .. Local Arrays ..
      DOUBLE PRECISION  C(17)
C     .. External Functions ..
      DOUBLE PRECISION  G02HKY
      EXTERNAL          G02HKY
C     .. External Subroutines ..
      EXTERNAL          C05AZF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
      XP = NVAR
      EP = 1.0D0/(1.0D0-EPS)
      A = 1.0D0
      B = 1.0D0
   20 FB = G02HKY(B,NVAR,TOL) - EP
      IF (FB.GE.0.0D0) THEN
         A = B
         B = B + 1.0D0
         GO TO 20
      END IF
   40 FA = G02HKY(A,NVAR,TOL) - EP
      IF (FA.LE.0.0D0) THEN
         B = A
         A = A/2.0D0
         GO TO 40
      END IF
      C(1) = FB
      IND = -1
      TOLX = TOL*XP
      IFAULT = 1
   60 CALL C05AZF(A,B,FA,TOLX,1,C,IND,IFAULT)
      IF (IND.NE.0) THEN
         FA = G02HKY(A,NVAR,TOL) - EP
         GO TO 60
      END IF
      A2 = MAX(XP-A,0.0D0)
      B2 = XP + A
      RETURN
      END
