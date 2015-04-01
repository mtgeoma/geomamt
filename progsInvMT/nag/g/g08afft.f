      SUBROUTINE G08AFF(X,LX,L,K,W1,H,P,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-336 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1118 (JUL 1993).
C     G08AFF PERFORMS THE KRUSKAL-WALLIS ONE-WAY ANALYSIS OF
C     VARIANCE BY RANKS ON K INDEPENDENT SAMPLES OF
C     POSSIBLY UNEQUAL SIZE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08AFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, P
      INTEGER           IFAIL, K, LX
C     .. Array Arguments ..
      DOUBLE PRECISION  W1(LX), X(LX)
      INTEGER           L(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIV, RJ, RS, T, XN, XNI
      INTEGER           I, I1, IDF, IFA, IFA2, J, L1, L2, LSUM
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08AEZ
C     .. Executable Statements ..
      IFA = 1
      IF (K.LT.2) GO TO 180
      IFA = 2
      LSUM = 0
      DO 20 I = 1, K
         IF (L(I).LE.0) GO TO 180
         LSUM = LSUM + L(I)
   20 CONTINUE
      IFA = 3
      IF (LSUM.NE.LX) GO TO 180
      DO 40 I = 2, LX
         IF (X(I).NE.X(1)) GO TO 60
   40 CONTINUE
      IFA = 4
      GO TO 180
   60 CONTINUE
C     RANK TOTAL SAMPLE (T = TIES CORRECTION TERM)
      CALL G08AEZ(X,W1,LX,1,T)
C     FIND RANK SUM FOR EACH SAMPLE
      RS = 0.0D0
      DO 140 I = 1, K
         L1 = 1
         IF (I.EQ.1) GO TO 100
         I1 = I - 1
         DO 80 J = 1, I1
            L1 = L1 + L(J)
   80    CONTINUE
  100    L2 = L1 + L(I) - 1
         RJ = 0.0D0
         DO 120 J = L1, L2
            RJ = RJ + W1(J)
  120    CONTINUE
         XNI = L(I)
         RJ = RJ*RJ/XNI
         RS = RS + RJ
  140 CONTINUE
C     COMPUTE TEST STATISTIC H
      XN = LX
      H = 12.0D0*RS/(XN*(XN+1.0D0)) - 3.0D0*(XN+1.0D0)
      IF (T.EQ.0.0D0) GO TO 160
      DIV = 1.0D0 - 12.0D0*T/(XN*XN*XN-XN)
      IFA = 4
      IF (DIV.LE.0.0D0) GO TO 180
      H = H/DIV
C     COMPUTE DEGREES OF FREEDOM AND SIGNIFICANCE
  160 IDF = K - 1
      IFA2 = 1
      P = G01ECF('U',H,DBLE(IDF),IFA2)
      IFAIL = 0
      GO TO 200
  180 IFAIL = P01ABF(IFAIL,IFA,SRNAME,0,P01REC)
  200 RETURN
      END
