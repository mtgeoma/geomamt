      SUBROUTINE G08AEF(X,IX,K,N,W1,W2,FR,P,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-335 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER 1117 (JUL 1993).
C     G08AEF PERFORMS THE FRIEDMAN TWO-WAY ANALYSIS OF VARIANCE
C     BY RANKS ON K RELATED SAMPLES OF SIZE N.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08AEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FR, P
      INTEGER           IFAIL, IX, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W1(K), W2(K), X(IX,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SSW, WMEAN, WW, XF
      INTEGER           I, IDF, IFA, IFA2, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08AEZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      IFA = 1
      IF (N.LT.1) GO TO 100
      IFA = 2
      IF (IX.LT.K) GO TO 100
      IFA = 3
      IF (K.LT.2) GO TO 100
C     RANK COLUMN 1 DIRECTLY INTO RANK SUM ARRAY W2
      CALL G08AEZ(X(1,1),W2,K,0,XF)
C     RANK COLUMNS 2 TO N INTO W1 THEN ADD TO W2
      IF (N.EQ.1) GO TO 60
      DO 40 J = 2, N
         CALL G08AEZ(X(1,J),W1,K,0,XF)
         DO 20 I = 1, K
            W2(I) = W2(I) + W1(I)
   20    CONTINUE
   40 CONTINUE
C     CALCULATE TOTAL SS ABOUT MEAN FOR ROW SUMS
   60 SSW = 0.0D0
      WMEAN = 0.5D0*DBLE(N*(K+1))
      DO 80 I = 1, K
         WW = W2(I) - WMEAN
         SSW = SSW + WW*WW
   80 CONTINUE
C     CALCULATE FRIEDMAN TEST VALUE,DEGREES OF FREEDOM,AND
C     SIGNIFICANCE
      FR = 12.0D0*SSW/DBLE(N*K*(K+1))
      IDF = K - 1
      IFA2 = 1
      P = G01ECF('U',FR,DBLE(IDF),IFA2)
      IF (P.GT.1.0D0) P = 1.0D0
      IFAIL = 0
      GO TO 120
  100 IFAIL = P01ABF(IFAIL,IFA,SRNAME,0,P01REC)
  120 RETURN
      END
