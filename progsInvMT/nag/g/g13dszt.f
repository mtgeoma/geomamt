      SUBROUTINE G13DSZ(IK,K,K2,M,N,R,P,Q,PARHLD,NPAR,R0,TEMP,WORK,SUM,
     *                  IDF,SIGLEV,IFAILA,IFAULT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15 REVISED. IER-926 (APR 1991).
C     MARK 16 REVISED. IER-1122 (JUL 1993).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGLEV, SUM
      INTEGER           IDF, IFAILA, IFAULT, IK, K, K2, M, N, NPAR, P, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  R(IK,IK,M), R0(K,K), TEMP(K2+1,K2), WORK(K2)
      LOGICAL           PARHLD(NPAR)
C     .. Local Scalars ..
      INTEGER           I, I2, IFAIL, J, J2, L
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      EXTERNAL          G01ECF
C     .. External Subroutines ..
      EXTERNAL          F01ADF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     This subroutine calculates the portmanteau statistic and its'
C     significance level
C
      IFAILA = 0
      IF (IFAULT.EQ.3) GO TO 280
C
C     first copy R0 onto TEMP
C
      DO 40 J = 1, K
         DO 20 I = 1, K
            TEMP(I,J) = R0(I,J)
   20    CONTINUE
   40 CONTINUE
C
C     invert RO
C
      IFAIL = 1
      CALL F01ADF(K,TEMP,K2+1,IFAIL)
      IF (IFAIL.GT.0) THEN
         IFAILA = 3
         GO TO 280
      END IF
C                -1
C     set R0 = R0
C
      DO 80 I = 1, K
         DO 60 J = 1, I
            R0(I,J) = TEMP(I+1,J)
            R0(J,I) = R0(I,J)
   60    CONTINUE
   80 CONTINUE
C               -1     -1
C     compute R0  ** R0   where ** denotes kronecker product
C
      DO 160 J = 1, K
         DO 140 I = 1, K
            SUM = R0(I,J)
            DO 120 J2 = 1, K
               DO 100 I2 = 1, K
                  TEMP((I-1)*K+I2,(J-1)*K+J2) = SUM*R0(I2,J2)
  100          CONTINUE
  120       CONTINUE
  140    CONTINUE
  160 CONTINUE
C
C     calculate Li-McLeod statistic
C
      SUM = 0.0D0
C
      DO 260 L = 1, M
C
C        copy vec(R(L)') onto WORK
C
C        Note that these R(L)'s are the transpose of those defined in
C        Li and McLeod's paper
C
         DO 200 I = 1, K
            DO 180 J = 1, K
               WORK((I-1)*K+J) = R(J,I,L)
  180       CONTINUE
  200    CONTINUE
C
         DO 240 J = 1, K2
            DO 220 I = 1, K2
               SUM = SUM + WORK(I)*TEMP(I,J)*WORK(J)
  220       CONTINUE
  240    CONTINUE
C
  260 CONTINUE
C
      SUM = SUM*DBLE(N)
C
C     modify Li-McLeod statistic
C
      SUM = SUM + (DBLE(K2*M*(M+1))/DBLE(2*N))
C
C     calculate degrees of freedom
C
  280 IDF = K2*(M-P-Q)
      DO 300 I = 1, NPAR
         IF (PARHLD(I)) IDF = IDF + 1
  300 CONTINUE
      IF (IFAILA.EQ.3 .OR. IFAULT.EQ.3) RETURN
      IFAIL = 1
      SIGLEV = G01ECF('U',SUM,DBLE(IDF),IFAIL)
      RETURN
C
      END
