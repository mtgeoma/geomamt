      DOUBLE PRECISION FUNCTION G05DFF(A,B)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER FROM THE CAUCHY DISTRIBUTION WITH
C     MEDIAN A AND SEMI-INTERQUARTILE RANGE B.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
C     .. Local Scalars ..
      DOUBLE PRECISION                 ONE, TWO, X, Y
C     .. External Functions ..
      DOUBLE PRECISION                 G05CAF
      EXTERNAL                         G05CAF
C     .. Data statements ..
      DATA                             ONE/1.0D0/, TWO/2.0D0/
C     .. Executable Statements ..
   20 X = G05CAF(X)*TWO - ONE
      Y = G05CAF(X)
      IF ((X*X+Y*Y).GT.ONE) GO TO 20
      G05DFF = A + B*X/Y
      RETURN
      END
