      SUBROUTINE G05EDF(N,P,R,NR,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SETS UP THE REFERENCE VECTOR FOR A BINOMIAL
C     DISTRIBUTION.
C     ADD AND TIMES CAN BE CHANGED IF A TRUNCATION PROBABILITY OF
C     1.0E-12 IS REGARDED AS UNSATISFACTORY.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IFAIL, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ADD, AV, CONST1, CONST2, CONST3, CONST4, HALF,
     *                  ONE, T, TIMES, TWOPI, U, V, W, X, Y, Z, ZERO
      INTEGER           I, IBASE, IBOT, IDIFF, IERR, ITOP, J, K
      LOGICAL           LSKEW
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05EXZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, DBLE, SQRT, INT
C     .. Data statements ..
      DATA              ADD/8.5D0/, TIMES/7.15D0/, ONE/1.0D0/,
     *                  ZERO/0.0D0/, HALF/0.5D0/,
     *                  TWOPI/6.283185307179586D0/,
     *                  CONST1/8.333333333333333D-2/,
     *                  CONST2/2.777777777777778D-3/,
     *                  CONST3/7.936507936507937D-4/,
     *                  CONST4/5.952380952380952D-4/
C     .. Executable Statements ..
      IERR = 1
      IF (N.LT.0) GO TO 240
      IERR = 2
      IF ((P.LT.ZERO) .OR. (P.GT.ONE)) GO TO 240
      AV = DBLE(N)*P
      Y = TIMES*SQRT(AV*(ONE-P))
      U = AV + Y + ADD
      V = AV - Y
      LSKEW = P .LE. HALF
      IF (LSKEW) GO TO 20
      U = U + ONE - ADD
      V = V + ONE - ADD
   20 ITOP = MIN(N,INT(U))
      IBOT = MAX(0,INT(V))
      IERR = 3
      IF ((ITOP-IBOT+4).GT.NR) GO TO 240
      IDIFF = ITOP - NR
      IBASE = IBOT - IDIFF
      W = DBLE(N+1)
      IF (( .NOT. LSKEW) .OR. (IBOT.GT.0)) GO TO 60
C     USE THE DIRECT METHOD (SKEW POSITIVE) IF NP .LT. 50(1-P).
      V = P/(ONE-P)
      X = ZERO
      Y = (ONE-P)**N
      Z = ZERO
      DO 40 I = IBASE, NR
         Z = Z + Y
         R(I) = Z
         X = X + ONE
         Y = Y*V*(W-X)/X
   40 CONTINUE
      GO TO 220
C     USE THE DIRECT METHOD (SKEW NEGATIVE) IF N(1-P) .LT. 50P.
   60 IF (LSKEW .OR. (ITOP.LT.N)) GO TO 120
      V = (ONE-P)/P
      X = ZERO
      Y = P**N
      DO 80 J = IBASE, NR
         I = NR + IBASE - J
         R(I) = Y
         X = X + ONE
         Y = Y*V*(W-X)/X
   80 CONTINUE
      Z = ZERO
      DO 100 I = IBASE, NR
         Z = Z + R(I)
         R(I) = Z
  100 CONTINUE
      GO TO 220
C     USE STIRLINGS FORMULA IF NEITHER.
  120 J = INT(AV)
      X = DBLE(J)
      Y = DBLE(N)
      Z = ONE/(Y*Y)
      V = (CONST1-(CONST2-CONST3*Z)*Z)/Y
      Z = ONE/(X*X)
      V = V - (CONST1-(CONST2-(CONST3-CONST4*Z)*Z)*Z)/X
      U = Y - X
      Z = ONE/(U*U)
      V = V - (CONST1-(CONST2-(CONST3-CONST4*Z)*Z)*Z)/U
      U = ZERO
      Z = ONE
      T = ONE
C     THIS IS EXP FOR SUITABLY SMALL ARGUMENTS.
      DO 140 I = 1, 5
         U = U + ONE
         Z = Z*V/U
         T = T + Z
  140 CONTINUE
      T = T*((P*Y/X)**J)*(((ONE-P)*Y/(Y-X))**(N-J))*SQRT(Y/(X*(Y-X)
     *    *TWOPI))
      J = J - IDIFF
      V = X
      Z = T
      U = (ONE-P)/P
      DO 160 K = IBASE, J
         I = J + IBASE - K
         R(I) = T
         T = T*U*V/(W-V)
         V = V - ONE
  160 CONTINUE
      Y = ZERO
      DO 180 I = IBASE, J
         Y = Y + R(I)
         R(I) = Y
  180 CONTINUE
      J = J + 1
      DO 200 I = J, NR
         X = X + ONE
         Z = Z*(W-X)/(U*X)
         Y = Y + Z
         R(I) = Y
  200 CONTINUE
C     FINISH OFF IN ALL CASES.
  220 CALL G05EXZ(IDIFF,IBASE,R,NR)
      IFAIL = 0
      RETURN
  240 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
