      SUBROUTINE G05ECF(T,R,NR,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SETS UP THE REFERENCE VECTOR FOR A POISSON DISTRIBUTION.
C     ADD AND TIMES CAN BE CHANGED IF A TRUNCATION PROBABILITY OF
C     1.0E-12 IS REGARDED AS UNSATISFACTORY.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05ECF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IFAIL, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ADD, CONST1, CONST2, CONST3, ONE, TIMES, TWOPI,
     *                  W, X, Y, Z, ZERO
      INTEGER           I, IBASE, IBOT, IDIFF, IERR, ITOP, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05EXZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, EXP, DBLE, SQRT, INT
C     .. Data statements ..
      DATA              ADD/8.5D0/, TIMES/7.15D0/, ZERO/0.0D0/,
     *                  ONE/1.0D0/, TWOPI/6.283185307179586D0/,
     *                  CONST1/8.333333333333333D-2/,
     *                  CONST2/2.777777777777778D-3/,
     *                  CONST3/7.936507936507937D-4/
C     .. Executable Statements ..
      IERR = 1
      IF (T.LT.ZERO) GO TO 140
      X = SQRT(T)*TIMES
      ITOP = INT(T+X+ADD)
      IBOT = MAX(0,INT(T-X))
      IERR = 2
      IF ((ITOP-IBOT+4).GT.NR) GO TO 140
      IDIFF = ITOP - NR
      IBASE = IBOT - IDIFF
      IF (IBOT.GT.0) GO TO 40
C     USE THE DIRECT METHOD IF T .LT. 50.
      X = ZERO
      Y = EXP(-T)
      Z = ZERO
      DO 20 I = IBASE, NR
         Z = Z + Y
         R(I) = Z
         X = X + ONE
         Y = Y*T/X
   20 CONTINUE
      GO TO 120
C     USE STIRLINGS FORMULA IF T .GT. 50.
   40 I = INT(T)
      X = DBLE(I)
      Z = ONE/(X*X)
      Y = ((T/X)**I)*EXP((X-T)-(CONST1-(CONST2-CONST3*Z)*Z)/X)
     *    /SQRT(TWOPI*X)
      J = I - IDIFF
      W = X
      Z = Y
      DO 60 K = IBASE, J
         I = J + IBASE - K
         R(I) = Y
         Y = Y*X/T
         X = X - ONE
   60 CONTINUE
      Y = ZERO
      DO 80 I = IBASE, J
         Y = Y + R(I)
         R(I) = Y
   80 CONTINUE
      J = J + 1
      DO 100 I = J, NR
         W = W + ONE
         Z = Z*T/W
         Y = Y + Z
         R(I) = Y
  100 CONTINUE
C     FINISH OFF IN EITHER CASE.
  120 CALL G05EXZ(IDIFF,IBASE,R,NR)
      IFAIL = 0
      RETURN
  140 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
