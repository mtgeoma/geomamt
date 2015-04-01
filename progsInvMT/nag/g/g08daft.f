      SUBROUTINE G08DAF(X,IX,K,N,RNK,W,P,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER 1119 (JUL 1993).
C
C     G08DAF CALCULATES KENDALLS COEFFICIENT OF CONCORDANCE
C     ON K INDEPENDENT RANKINGS OF N OBJECTS OR INDIVIDUALS.
C
C     REFERENCE - S SIEGEL - NONPARAMETRIC STATISTICS FOR
C                            THE BEHAVIORAL SCIENCES ( P 229 )
C     AUTHOR    - J. LLOYD-JONES (U.M.R.C.C.)
C
C     PARAMETERS -
C     X     - MATRIX OF RANKS
C     IX    - FIRST DIMENSION OF X
C     K     - NUMBER OF COMPARISONS
C     N     - NUMBER OF OBJECTS BEING COMPARED
C     RNK   - RANKED DATA
C     W     - COEFFICIENT OF CONCORDANCE
C     P     - SIGNIFICANCE LEVEL OUTPUT BY G01ECF
C     IFAIL - ERROR PARAMETER
C
C     LOCAL VARIABLES -
C     RMEAN - OVERALL AVERAGE RANK
C     IDF   - DEGREES OF FREEDOM(=N-1) INPUT TO G01ECF
C     CHI   - CHI SQUARE VALUE INPUT TO G01ECF
C
C     CHECK PARAMETERS
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08DAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, W
      INTEGER           IFAIL, IX, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RNK(IX,N), X(IX,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CHI, FK, FN, R, RMEAN, T
      INTEGER           I, IDF, IERROR, IFA, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08DAZ
C     .. Executable Statements ..
      IF (N.LT.2) GO TO 60
      IF (IX.LT.K) GO TO 80
      IF (K.LT.2) GO TO 100
C     ********************************************************
C     CALL AUXILIARY G08DAZ TO RANK VALUES AND
C     CALCULATE T , THE CORRECTION FACTOR FOR TIED VALUES
C     ********************************************************
      CHI = 0.D0
      CALL G08DAZ(X,IX,K,N,RNK,T)
C     ************
C     CALCULATE W
C     ************
      W = 0.D0
      FK = K
      FN = N
      RMEAN = 0.5D0*FK*(FN+1.0D0)
      DO 40 J = 1, N
         R = RNK(1,J)
         DO 20 I = 2, K
            R = R + RNK(I,J)
   20    CONTINUE
         W = W + (R-RMEAN)*(R-RMEAN)
   40 CONTINUE
      W = W/(FK*FK*FN*(FN*FN-1.0D0)/12.0D0-FK*T)
C     **************
C     SIGNIFICANCE
C     **************
      CHI = W*FK*(FN-1.0D0)
      IDF = N - 1
      IFA = 1
      P = G01ECF('U',CHI,DBLE(IDF),IFA)
      IFAIL = 0
      GO TO 140
   60 IERROR = 1
      GO TO 120
   80 IERROR = 2
      GO TO 120
  100 IERROR = 3
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
  140 RETURN
      END
