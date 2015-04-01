      SUBROUTINE G05EFF(L,M,N,R,NR,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SETS UP THE REFERENCE VECTOR FOR A HYPERGEOMETRIC
C     DISTRIBUTION.
C     ADD AND TIMES CAN BE CHANGED IF A TRUNCATION PROBABILITY OF
C     1.0E-12 IS REGARDED AS UNSATISFACTORY.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, L, M, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ADD, AK, AL, AM, AN, CONST1, CONST2, CONST3,
     *                  CONST4, CONST5, ONE, T, TIMES, TWOPI, U, V, W,
     *                  X, Y, Z, ZERO
      INTEGER           I, IBASE, IBOT, IDIFF, IERR, IMAX, IMIN, ITOP,
     *                  J, J1, K, K1, K2
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
      DATA              ADD/8.5D0/, TIMES/7.15D0/, ZERO/0.0D0/,
     *                  ONE/1.0D0/, TWOPI/6.283185307179586D0/,
     *                  CONST1/8.333333333333333D-2/,
     *                  CONST2/2.777777777777778D-3/,
     *                  CONST3/7.936507936507937D-4/,
     *                  CONST4/5.952380952380952D-4/,
     *                  CONST5/8.417508417508418D-4/
C     .. Executable Statements ..
      IERR = 1
      IF (N.LT.0) GO TO 480
      IERR = 2
      IF ((L.LT.0) .OR. (L.GT.N)) GO TO 480
      IERR = 3
      IF ((M.LT.0) .OR. (M.GT.N)) GO TO 480
      IMIN = MAX((L+M-N),0)
      IMAX = MIN(L,M)
      LSKEW = (((L+L).LE.N) .AND. ((M+M).LE.N)) .OR. (((L+L).GE.N)
     *        .AND. ((M+M).GE.N))
      AL = DBLE(L)
      AM = DBLE(M)
      AN = DBLE(N)
      IF (N.GT.1) GO TO 20
      AK = AL*AM
      X = ZERO
      GO TO 40
   20 AK = AL*AM/AN
      X = TIMES*SQRT(AK*(AN-AL)*(AN-AM)/(AN*(AN-ONE)))
   40 Y = AK - X
      Z = AK + X + ADD
      IF (LSKEW) GO TO 60
      Y = Y + ONE - ADD
      Z = Z + ONE - ADD
   60 IBOT = MAX(IMIN,INT(Y))
      ITOP = MIN(IMAX,INT(Z))
      IERR = 4
      IF ((ITOP-IBOT+4).GT.NR) GO TO 480
      IDIFF = ITOP - NR
      IBASE = IBOT - IDIFF
      K = INT(AK)
      AK = DBLE(K)
      IF (( .NOT. LSKEW) .OR. (IBOT.GT.IMIN)) GO TO 200
C     USE THE DIRECT METHOD (SKEW POSITIVE) IF
C     LM/N .LT. 50(N-L)(N-M)/(N*N) OR (N-L)(N-M)/N .LT. 50LM/(N*N).
      IF (IMIN.GT.0) GO TO 80
      K1 = L
      K2 = M
      GO TO 100
   80 K1 = N - L
      K2 = N - M
  100 IF (K1.LE.K2) GO TO 120
      I = K1
      K1 = K2
      K2 = I
  120 Y = ONE
      U = DBLE(K1)
      V = DBLE(K2)
      IF (K1.LE.0) GO TO 160
      W = AN - U
      DO 140 I = 1, K1
         W = W + ONE
         Y = Y*(W-V)/W
  140 CONTINUE
  160 Z = ZERO
      W = ONE
      X = AN - U - V
      U = U + ONE
      V = V + ONE
      DO 180 I = IBASE, NR
         Z = Z + Y
         R(I) = Z
         Y = Y*(U-W)*(V-W)/(W*(X+W))
         W = W + ONE
  180 CONTINUE
      GO TO 460
  200 IF (LSKEW .OR. (ITOP.LT.IMAX)) GO TO 360
C     USE THE DIRECT METHOD (SKEW NEGATIVE) IF
C     L(N-M)/N .LT. 50(N-L)M/(N*N) OR (N-L)M/N .LT. 50L(N-M)/(N*N).
      IF (IMAX.LT.L) GO TO 220
      K1 = L
      K2 = N - M
      GO TO 240
  220 K1 = N - L
      K2 = M
  240 IF (K1.LE.K2) GO TO 260
      I = K1
      K1 = K2
      K2 = I
  260 Y = ONE
      U = DBLE(K1)
      V = DBLE(K2)
      IF (K1.LE.0) GO TO 300
      W = AN - U
      DO 280 I = 1, K1
         W = W + ONE
         Y = Y*(W-V)/W
  280 CONTINUE
  300 W = ONE
      X = AN - U - V
      U = U + ONE
      V = V + ONE
      DO 320 I = IBASE, NR
         J = NR + IBASE - I
         R(J) = Y
         Y = Y*(U-W)*(V-W)/(W*(X+W))
         W = W + ONE
  320 CONTINUE
      Z = ZERO
      DO 340 I = IBASE, NR
         Z = Z + R(I)
         R(I) = Z
  340 CONTINUE
      GO TO 460
C     USE STIRLINGS FORMULA WHEN NONE IS TRUE.
  360 U = ONE/(AL*AL)
      W = (CONST1-(CONST2-(CONST3-CONST4*U)*U)*U)/AL
      U = ONE/(AM*AM)
      W = W + (CONST1-(CONST2-(CONST3-CONST4*U)*U)*U)/AM
      V = AN - AL
      U = ONE/(V*V)
      W = W + (CONST1-(CONST2-(CONST3-CONST4*U)*U)*U)/V
      V = AN - AM
      U = ONE/(V*V)
      W = W + (CONST1-(CONST2-(CONST3-CONST4*U)*U)*U)/V
      U = ONE/(AK*AK)
      W = W - (CONST1-(CONST2-(CONST3-(CONST4-CONST5*U)*U)*U)*U)/AK
      V = AL - AK
      U = ONE/(V*V)
      W = W - (CONST1-(CONST2-(CONST3-(CONST4-CONST5*U)*U)*U)*U)/V
      V = AM - AK
      U = ONE/(V*V)
      W = W - (CONST1-(CONST2-(CONST3-(CONST4-CONST5*U)*U)*U)*U)/V
      V = AN - AL - AM + AK
      U = ONE/(V*V)
      W = W - (CONST1-(CONST2-(CONST3-(CONST4-CONST5*U)*U)*U)*U)/V
      U = ONE/(AN*AN)
      W = W - (CONST1-(CONST2-CONST3*U)*U)/AN
      U = ZERO
      V = ONE
      X = ONE
C     THIS IS EXP FOR SUITABLY SMALL ARGUMENTS.
      DO 380 I = 1, 6
         U = U + ONE
         V = V*W/U
         X = X + V
  380 CONTINUE
      X = X*SQRT(AL*AM*(AN-AL)*(AN-AM)/(AK*AN*(AL-AK)*(AM-AK)
     *    *(AN-AL-AM+AK)*TWOPI))*(AL*AM/(AK*AN))**K*(AL*(AN-AM)/((AL-AK)
     *    *AN))**(L-K)*(AM*(AN-AL)/((AM-AK)*AN))**(M-K)*((AN-AL)*(AN-AM)
     *    /((AN-AL-AM+AK)*AN))**(N-L-M+K)
      J = K - IDIFF
      U = AK
      T = X
      W = AM + ONE
      Y = AL + ONE
      V = AN - AM - AL
      DO 400 I = IBASE, J
         J1 = J + IBASE - I
         R(J1) = X
         X = X*U*(V+U)/((W-U)*(Y-U))
         U = U - ONE
  400 CONTINUE
      Z = ZERO
      DO 420 I = IBASE, J
         Z = Z + R(I)
         R(I) = Z
  420 CONTINUE
      J = J + 1
      U = AK
      DO 440 I = J, NR
         U = U + ONE
         T = T*(W-U)*(Y-U)/(U*(V+U))
         Z = Z + T
         R(I) = Z
  440 CONTINUE
C     FINISH OFF IN ALL CASES.
  460 CALL G05EXZ(IDIFF,IBASE,R,NR)
      IFAIL = 0
      RETURN
  480 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
