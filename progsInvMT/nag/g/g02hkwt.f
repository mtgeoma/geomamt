      SUBROUTINE G02HKW(A2,B2,NVAR,T2)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C      CALCULATES THE CORRECTION TERM TAU**2
C      FOR ROBUST COVARIANCE ESTIMATION USING
C      HUBER'S WEIGHT FUNCTION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A2, B2, T2
      INTEGER           NVAR
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, FB, TOL, XP
      INTEGER           IFAULT, IND
C     .. Local Arrays ..
      DOUBLE PRECISION  C(17)
C     .. External Functions ..
      DOUBLE PRECISION  G02HKV
      EXTERNAL          G02HKV
C     .. External Subroutines ..
      EXTERNAL          C05AZF
C     .. Executable Statements ..
      TOL = 0.5D-5
      XP = NVAR
      A = 0.0D0
      B = 1.0D0
   20 FB = G02HKV(B,XP,A2,B2) - XP
      IF (FB.LE.0.0D0) THEN
         A = B
         B = B + 1.0D0
         GO TO 20
      END IF
      C(1) = A2 - XP
      IND = -1
      IFAULT = 1
   40 CALL C05AZF(B,A,FB,TOL,2,C,IND,IFAULT)
      IF (IND.NE.0) THEN
         FB = G02HKV(B,XP,A2,B2) - XP
         GO TO 40
      END IF
      T2 = B
      RETURN
      END
