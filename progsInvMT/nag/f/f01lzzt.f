      DOUBLE PRECISION FUNCTION F01LZZ(A,B,SMALL,BIG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (TANGNT)
C
C     F01LZZ RETURNS THE VALUE
C
C     F01LZZ = B/A .
C
C     SMALL AND BIG MUST BE SUCH THAT
C
C     SMALL = X02AMF     AND     BIG = 1.0/SMALL ,
C
C     WHERE X02AMF IS THE SMALL NUMBER RETURNED FROM ROUTINE
C     X02AMF.
C
C     IF B/A IS LESS THAN SMALL THEN F01LZZ IS RETURNED AS
C     ZERO AND IF B/A IS GREATER THAN BIG THEN F01LZZ IS
C     RETURNED AS SIGN(BIG,B).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, BIG, SMALL
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSA, ABSB, X
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SIGN
C     .. Executable Statements ..
      F01LZZ = 0.0D0
      IF (B.EQ.0.0D0) RETURN
C
      ABSA = ABS(A)
      ABSB = ABS(B)
      X = 0.0D0
      IF (ABSA.GE.1.0D0) X = ABSA*SMALL
C
      IF (ABSB.LT.X) RETURN
C
      X = 0.0D0
      IF (ABSB.GE.1.0D0) X = ABSB*SMALL
C
      IF (ABSA.LE.X) GO TO 20
C
      F01LZZ = B/A
      RETURN
C
   20 F01LZZ = SIGN(BIG,B)
      RETURN
      END