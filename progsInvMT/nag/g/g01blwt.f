      SUBROUTINE G01BLW(K,N,L,M,LPK)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes log(P(K)), where
C
C     P(K) = (L!M!(N-L)!(N-M)!)/(K!(L-K)!(M-K)!(N-L-M+K)!N!)
C
C     Intended ranges of the input arguments K, N, L, M:
C
C     2 .LE. N
C     0 .LT. L .LT. N
C     0 .LT. M .LT. N
C     MAX(0,L+M-N) .LE. K .LE. MIN(L,M)
C
C     G01BJU returns log(1+X) for X.ge.0
C     G01BJV returns log(gamma(X+1)*exp(X)/X**X) for X.ge.0
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  LPK
      INTEGER           K, L, M, N
C     .. Local Scalars ..
      DOUBLE PRECISION  D, D2, RK, RL, RM, RN, T1, T2, T3, T4
      INTEGER           IFAIL
C     .. Local Arrays ..
      DOUBLE PRECISION  AV(2), BV(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01BJU, G01BJV
      EXTERNAL          G01BJU, G01BJV
C     .. External Subroutines ..
      EXTERNAL          X03AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      RN = DBLE(N)
      RL = DBLE(L)
      RM = DBLE(M)
      RK = DBLE(K)
C
C     Compute K*N - L*M, using additional precision in implementations
C     of low precision. In the call of X03AAF, the argument SW may be
C     changed to .FALSE. in double precision implementations.
C
      AV(1) = RK
      AV(2) = -RL
      BV(1) = RN
      BV(2) = RM
      IFAIL = 0
      CALL X03AAF(AV,2,BV,2,2,1,1,0.0D0,0.0D0,D,D2,.FALSE.,IFAIL)
C
      IF (K.EQ.0) THEN
         T1 = 0.0D0
      ELSE IF (D.GE.0.0D0) THEN
         T1 = -RK*G01BJU(D/(RL*RM))
      ELSE
         T1 = RK*G01BJU(-D/(RK*RN))
      END IF
C
      IF (L-K.EQ.0) THEN
         T2 = 0.0D0
      ELSE IF (D.GE.0.0D0) THEN
         T2 = (RL-RK)*G01BJU(D/(RN*(RL-RK)))
      ELSE
         T2 = -(RL-RK)*G01BJU(-D/(RL*(RN-RM)))
      END IF
C
      IF (M-K.EQ.0) THEN
         T3 = 0.0D0
      ELSE IF (D.GE.0.0D0) THEN
         T3 = (RM-RK)*G01BJU(D/(RN*(RM-RK)))
      ELSE
         T3 = -(RM-RK)*G01BJU(-D/(RM*(RN-RL)))
      END IF
C
      IF (N-L-M+K.EQ.0) THEN
         T4 = 0.0D0
      ELSE IF (D.GE.0.0D0) THEN
         T4 = -(RN-RL-RM+RK)*G01BJU(D/((RN-RL)*(RN-RM)))
      ELSE
         T4 = (RN-RL-RM+RK)*G01BJU(-D/(RN*(RN-RL-RM+RK)))
      END IF
C
      LPK = G01BJV(RL) + G01BJV(RM) + G01BJV(RN-RL) + G01BJV(RN-RM)
      LPK = LPK - G01BJV(RK) - G01BJV(RL-RK) - G01BJV(RM-RK) -
     *      G01BJV(RN-RL-RM+RK) - G01BJV(RN)
      LPK = LPK + T1 + T2 + T3 + T4
      RETURN
      END
