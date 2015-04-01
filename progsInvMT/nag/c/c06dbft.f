      DOUBLE PRECISION FUNCTION C06DBF(X,C,N,S)
C     MARK 6 RELEASE. NAG COPYRIGHT 1977.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS FUNCTION EVALUATES A CHEBYSHEV SERIES OF ONE OF THREE
C     FORMS ACCORDING TO THE VALUE OF S
C     S=1- 0.5*C(1)+SUM FROM J=2 TO N OF C(J)*T(J-1)(X),
C     S=2- 0.5*C(1)+SUM FROM J=2 TO N OF C(J)*T(2J-2)(X),
C     S=3- SUM FROM J=1 TO N OF C(J)*T(2J-1)(X),
C     WHERE X LIES IN THE RANGE -1.0 .LE. X .LE. 1.0
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          N, S
C     .. Array Arguments ..
      DOUBLE PRECISION                 C(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, B, D, Y, Y2
      INTEGER                          R
C     .. Executable Statements ..
      IF (S.EQ.1) GO TO 20
      Y = 2.0D0*X*X - 1.0D0
      GO TO 40
   20 Y = X
   40 Y2 = 2.0D0*Y
      A = 0.0D0
      B = 0.0D0
      R = N
      D = C(R)
      IF (R.EQ.1) GO TO 80
   60 R = R - 1
      A = B
      B = D
      D = C(R) - A + Y2*B
      IF (R.GT.1) GO TO 60
   80 IF (S.EQ.3) GO TO 100
      C06DBF = 0.5D0*(D-A)
      GO TO 120
  100 C06DBF = X*(D-B)
  120 RETURN
      END
