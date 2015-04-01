      SUBROUTINE G05EEF(N,P,R,NR,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SETS UP THE REFERENCE VECTOR FOR A NEGATIVE BINOMIAL
C     DISTRIBUTION.
C     ADD, TIMES AND TIMES2 CAN BE CHANGED IF A TRUNCATION
C     PROBABILITY
C     OF 1.0E-12 IS REGARDED AS UNSATISFACTORY.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IFAIL, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ADD, AV, CONST1, CONST2, CONST3, ONE, T, TIMES,
     *                  TIMES2, TWOPI, U, V, W, X, Y, Z, ZERO
      INTEGER           I, IBASE, IBOT, IDIFF, IERR, ITOP, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05EXZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT, INT
C     .. Data statements ..
      DATA              ADD/8.5D0/, TIMES/7.15D0/, TIMES2/20.15D0/,
     *                  ONE/1.0D0/, ZERO/0.0D0/,
     *                  TWOPI/6.283185307179586D0/,
     *                  CONST1/8.333333333333333D-2/,
     *                  CONST2/2.777777777777778D-3/,
     *                  CONST3/7.936507936507937D-4/
C     .. Executable Statements ..
      IERR = 1
      IF (N.LT.0) GO TO 160
      IERR = 2
      IF ((P.LT.ZERO) .OR. (P.GE.ONE)) GO TO 160
      AV = DBLE(N)*P/(ONE-P)
      X = TIMES*SQRT(AV/(ONE-P))
      ITOP = INT(AV+X+ADD+TIMES2*P/(ONE-P))
      IBOT = MAX(0,INT(AV-X))
      IERR = 3
      IF ((ITOP-IBOT+4).GT.NR) GO TO 160
      IDIFF = ITOP - NR
      IBASE = IBOT - IDIFF
      W = DBLE(N-1)
      IF (IBOT.GT.0) GO TO 40
C     USE THE DIRECT METHOD IF NP .LT. 50.
      X = ZERO
      Y = (ONE-P)**N
      Z = ZERO
      DO 20 I = IBASE, NR
         Z = Z + Y
         R(I) = Z
         X = X + ONE
         Y = Y*P*(W+X)/X
   20 CONTINUE
      GO TO 140
C     USE STIRLINGS FORMULA IF NP .GT. 50.
   40 J = INT(AV)
      X = DBLE(J)
      Y = X + W
      Z = ONE/(Y*Y)
      V = (CONST1-(CONST2-CONST3*Z)*Z)/Y
      Z = ONE/(X*X)
      V = V - (CONST1-(CONST2-CONST3*Z)*Z)/X
      Z = ONE/(W*W)
      V = V - (CONST1-(CONST2-CONST3*Z)*Z)/W
      U = ZERO
      Z = ONE
      T = ONE
C     THIS IS EXP FOR SUITABLY SMALL ARGUMENTS.
      DO 60 I = 1, 4
         U = U + ONE
         Z = Z*V/U
         T = T + Z
   60 CONTINUE
      T = T*(ONE-P)*(((ONE-P)*(X+W)/W)**(N-1))*((P*(X+W)/X)**J)
     *    *SQRT((X+W)/(X*W*TWOPI))
      J = J - IDIFF
      V = X
      Z = T
      DO 80 K = IBASE, J
         I = J + IBASE - K
         R(I) = T
         T = T*V/((W+V)*P)
         V = V - ONE
   80 CONTINUE
      Y = ZERO
      DO 100 I = IBASE, J
         Y = Y + R(I)
         R(I) = Y
  100 CONTINUE
      J = J + 1
      DO 120 I = J, NR
         X = X + ONE
         Z = Z*P*(W+X)/X
         Y = Y + Z
         R(I) = Y
  120 CONTINUE
C     FINISH OFF IN EITHER CASE.
  140 CALL G05EXZ(IDIFF,IBASE,R,NR)
      IFAIL = 0
      RETURN
  160 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
