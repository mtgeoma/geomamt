      SUBROUTINE G08RAZ(ZIN,ETA,VAPVEC,N,N1,IDIST)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATES EXPECTED VALUES AND VARIANCE-COVARIANCE MATRIX
C     OF PARTICULAR FUNCTION OF ORDER STATISTICS FOR THE NORMAL,
C     LOGISTIC, DOUBLE EXPONENTIAL AND EXTREME VALUE DISTRIBUTIONS.
C
C     PETTITT A.N.P. INFERENCE FOR THE LINEAR MODEL USING A LIKELIHOOD
C                    BASED ON RANKS.
C                    JRSS,B,44 PP 234-243
C
C     ARGUMENTS :
C                 ZIN - EXPECTED VALUE OF FUNCTION OF ORDER STATISTIC
C                 ETA - EXPECTED VALUE OF DERIVATIVE OF FUNCTION OF
C                       ORDER STATISTIC.
C              VAPVEC - VARIANCE-COVARIANCE MATRIX OF FUNCTION OF ORDER
C                       STATISTIC.
C                   N - SAMPLE SIZE
C                  N1 - DIMENSION OF VECTOR NEEDED TO STORE UPPER
C                       TRIANGLE OF SQUARE SYMMETRIC MATRIX.
C                       ( N1 .GE. N*(N+1)/2 ).
C               IDIST - SPECIFIES ERROR DISTRIBUTION TO BE USED
C                       1 - NORMAL
C                       2 - LOGISTIC
C                       3 - DOUBLE EXPONENTIAL
C                       4 - EXTREME VALUE
C
C     .. Scalar Arguments ..
      INTEGER           IDIST, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), VAPVEC(N1), ZIN(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIV, FOUR, ONE, SUM, SUM1, SUM2, SUMM2, TERM,
     *                  TWO, V11, XL2, XNL2, ZERO
      INTEGER           I, IFAIL, J, M, N2
C     .. External Functions ..
      DOUBLE PRECISION  G01DCV
      INTEGER           G01DCU
      EXTERNAL          G01DCV, G01DCU
C     .. External Subroutines ..
      EXTERNAL          G01DBF, G01DCF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP
C     .. Data statements ..
      DATA              XL2/0.6931471805D0/, ZERO/0.0D0/, ONE/1.0D0/,
     *                  TWO/2.0D0/, FOUR/4.0D0/
C     .. Executable Statements ..
C
C     GO TO APPROPRIATE ERROR DISTRIBUTION
C
      GO TO (20,60,120,180) IDIST
C
C     FOR THE NORMAL DISTRIBUTION
C
   20 IFAIL = 0
      CALL G01DBF(N,ZIN,IFAIL)
      SUMM2 = ZERO
      DO 40 I = 1, N
         SUMM2 = SUMM2 + ZIN(I)*ZIN(I)
         ETA(I) = ONE
   40 CONTINUE
      V11 = G01DCV(N)
      IFAIL = 0
      CALL G01DCF(N,ZIN(N),ZIN(N-1),SUMM2,VAPVEC,IFAIL)
      RETURN
C
C     FOR THE LOGISTIC DISTRIBUTION
C
   60 DO 100 I = 1, N
         DO 80 J = 1, I
            M = G01DCU(I,J)
            VAPVEC(M) = FOUR*J*(N+ONE-I)/((N+1)**2*(N+2))
   80    CONTINUE
         ZIN(I) = TWO*I/(N+ONE) - ONE
         M = G01DCU(I,I)
         ETA(I) = (N+ONE)*VAPVEC(M)/TWO
  100 CONTINUE
      RETURN
C
C     FOR THE DOUBLE EXPONENTIAL DISTRIBUTION
C
  120 SUM1 = ZERO
      SUM2 = ZERO
      DO 160 J = 1, N
         SUM2 = SUM2 + ONE/((N-J+ONE)*(N-J+ONE))
         DO 140 I = J, N
            M = G01DCU(I,J)
            VAPVEC(M) = SUM2
  140    CONTINUE
         SUM1 = SUM1 + ONE/(N-J+ONE)
         ZIN(J) = SUM1 - ONE
         ETA(J) = SUM1
  160 CONTINUE
      RETURN
C
C     FOR THE EXTREME VALUE DISTRIBUTION
C
  180 XNL2 = N*XL2
      SUM = ZERO
      TERM = ONE
      N2 = N/2
      IF (2*N2.NE.N) ZIN(N2+1) = ZERO
      DO 200 I = 1, N2
         DIV = I
         SUM = SUM + TERM
         ZIN(N-I+1) = ONE - TWO*EXP(-XNL2)*SUM
         ETA(N-I+1) = TWO*(N+1-I)*EXP(-XNL2)*TERM
         ZIN(I) = -ZIN(N-I+1)
         ETA(I) = ETA(N-I+1)
         TERM = TERM*(N-I+ONE)/DIV
  200 CONTINUE
      IF (2*N2.EQ.N) GO TO 220
      ETA(N2+1) = TWO*(N2+ONE)*EXP(-XNL2)*TERM
  220 DO 260 J = 1, N
         DO 240 I = 1, J
            M = G01DCU(I,J)
            VAPVEC(M) = ONE - ZIN(J) + ZIN(I) - ZIN(I)*ZIN(J)
  240    CONTINUE
  260 CONTINUE
      RETURN
      END
