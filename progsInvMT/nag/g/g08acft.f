      SUBROUTINE G08ACF(X,N,N1,W1,I1,I2,P,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1116 (JUL 1993).
C     G08ACF PERFORMS THE MEDIAN TEST ON TWO INDEPENDENT
C     SAMPLES OF POSSIBLY UNEQUAL SIZE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08ACF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           I1, I2, IFAIL, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  W1(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CHIS, P1, P2, RM, XF
      INTEGER           I, IFA, IFA2, INOB, IPRED, M3, MM, N1P1, N2, N3,
     *                  NDF, NN, NPOS, NUM
C     .. Local Arrays ..
      DOUBLE PRECISION  PRED(3,3), PROB(21)
      INTEGER           NOBS(3,3)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01AFF, G08AEZ
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IFA = 1
      IF (N.LT.2) GO TO 140
      IFA = 2
      IF (N1.GE.N .OR. N1.LT.1) GO TO 140
      N2 = N - N1
      CALL G08AEZ(X,W1,N,0,XF)
C     RM = LARGEST RANK BELOW MEDIAN
      RM = (N+1)/2
C     I1 = NUMBER IN SAMPLE 1 WHICH ARE BELOW MEDIAN
C     I2 = NUMBER IN SAMPLE 2 WHICH ARE BELOW MEDIAN
      I1 = 0
      DO 20 I = 1, N1
         IF (W1(I).LT.RM) I1 = I1 + 1
   20 CONTINUE
      I2 = 0
      N1P1 = N1 + 1
      DO 40 I = N1P1, N
         IF (W1(I).LT.RM) I2 = I2 + 1
   40 CONTINUE
C     SET UP PARAMETERS FOR G01AFF
      INOB = 3
      IPRED = 3
      M3 = 3
      N3 = 3
      NOBS(1,1) = I1
      NOBS(1,2) = I2
      NOBS(2,1) = N1 - I1
      NOBS(2,2) = N2 - I2
      NUM = 0
      IFA2 = 1
      CALL G01AFF(INOB,IPRED,M3,N3,NOBS,NUM,PRED,CHIS,PROB,NPOS,NDF,MM,
     *            NN,IFA2)
      IF (NUM.GT.0) GO TO 60
C     COMPUTE CHI-SQUARE PROBABILITY
      IFA2 = 1
      P = G01ECF('U',CHIS,DBLE(NDF),IFA2)
      GO TO 120
C     COMPUTE FISHERS EXACT TEST PROBABILITY
   60 P1 = 0.0D0
      P2 = 0.0D0
      DO 80 I = 1, NPOS
         P1 = P1 + PROB(I)
   80 CONTINUE
      DO 100 I = NPOS, NUM
         P2 = P2 + PROB(I)
  100 CONTINUE
      P = MIN(P1,P2)
  120 IFAIL = 0
      GO TO 160
  140 IFAIL = P01ABF(IFAIL,IFA,SRNAME,0,P01REC)
  160 RETURN
      END
