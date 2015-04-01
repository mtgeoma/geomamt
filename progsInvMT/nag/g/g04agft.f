      SUBROUTINE G04AGF(Y,N,K,LSUB,NOBS,L,NGP,GBAR,SGBAR,GM,SS,IDF,F,FP,
     *                  IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1115 (JUL 1993).
C     G04AGF PERFORMS AN ANALYSIS OF VARIANCE FOR A
C     TWO-WAY HIERARCHICAL CLASSIFICATION WITH
C     SUBGROUPS OF POSSIBLY UNEQUAL SIZE , AND ALSO
C     COMPUTES THE TREATMENT GROUP AND
C     SUBGROUP MEANS.
C     CHECK PARAMETERS
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G04AGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GM
      INTEGER           IFAIL, K, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  F(2), FP(2), GBAR(K), SGBAR(L), SS(4), Y(N)
      INTEGER           IDF(4), LSUB(K), NGP(K), NOBS(L)
C     .. Local Scalars ..
      DOUBLE PRECISION  FDEN, FF, S1, S2, S4, YDD, YID, YIJ, YY, Z
      INTEGER           I, IFA, IFAULT, J, LI, LSUM, M, NGPI, NHI, NIJ,
     *                  NLO, NSUB, NSUM
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01EDF
      INTEGER           P01ABF
      EXTERNAL          G01EDF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      IFAULT = 1
      IF (K.LE.1) GO TO 260
      IFAULT = 2
      LSUM = 0
      DO 20 I = 1, K
         IF (LSUB(I).LE.0) GO TO 260
         LSUM = LSUM + LSUB(I)
   20 CONTINUE
      IFAULT = 3
      IF (LSUM.NE.L) GO TO 260
      IFAULT = 4
      NSUM = 0
      DO 40 I = 1, L
         IF (NOBS(I).LE.0) GO TO 260
         NSUM = NSUM + NOBS(I)
   40 CONTINUE
      IFAULT = 5
      IF (NSUM.NE.N) GO TO 260
      IFAULT = 6
C     COMPUTE   GRAND MEAN      GM
C     TOTAL SS ABOUT MEAN    S4
C     GP I SGP J MEANS
C     GROUP I MEANS
      NLO = 1
      NSUB = 0
      YDD = 0.0D0
      S4 = 0.0D0
      DO 100 I = 1, K
         NGPI = 0
         YID = 0.0D0
         LI = LSUB(I)
         DO 80 J = 1, LI
            YIJ = 0.0D0
            NSUB = NSUB + 1
            NIJ = NOBS(NSUB)
            NGPI = NGPI + NIJ
            NHI = NLO + NIJ - 1
            DO 60 M = NLO, NHI
               YY = Y(M)
               YDD = YDD + YY
               YID = YID + YY
               YIJ = YIJ + YY
   60       CONTINUE
            SGBAR(NSUB) = YIJ/DBLE(NIJ)
            NLO = NLO + NIJ
   80    CONTINUE
         NGP(I) = NGPI
         GBAR(I) = YID/DBLE(NGPI)
  100 CONTINUE
      GM = YDD/DBLE(N)
      DO 120 I = 1, N
         Z = Y(I) - GM
         S4 = S4 + Z*Z
  120 CONTINUE
C     ALL MEANS COMPUTED - EXIT IF TOTAL SS ABOUT MEAN IS ZERO
      DO 140 I = 1, 4
         SS(I) = 0.0D0
  140 CONTINUE
      IF (S4.LE.0.0D0) GO TO 260
      IFAULT = 7
C     COMPUTE    GROUP SS     S1 = SS(1)
C     SUBGROUP SS  S2 = SS(2)
      S1 = 0.0D0
      S2 = 0.0D0
      NSUB = 0
      DO 180 I = 1, K
         Z = GBAR(I) - GM
         S1 = S1 + Z*Z*DBLE(NGP(I))
         LI = LSUB(I)
         DO 160 J = 1, LI
            NSUB = NSUB + 1
            Z = SGBAR(NSUB) - GBAR(I)
            S2 = S2 + Z*Z*DBLE(NOBS(NSUB))
  160    CONTINUE
  180 CONTINUE
      SS(1) = S1
      SS(2) = S2
      SS(3) = S4 - S2 - S1
      SS(4) = S4
C     ASSIGN DEGREES OF FREEDOM
      IDF(1) = K - 1
      IDF(2) = L - K
      IDF(3) = N - L
      IDF(4) = N - 1
      IF (SS(3).LE.0.0D0) GO TO 260
C     COMPUTE F RATIOS AND SIGNIFICANCES
      FDEN = SS(3)/DBLE(IDF(3))
      DO 240 I = 1, 2
         FF = SS(I)/DBLE(IDF(I))
         F(I) = 9999.0D0
         FP(I) = 0.0D0
         IF (FDEN.GE.1.0D0) GO TO 200
         IF (FF.GE.FDEN*9999.0D0) GO TO 240
         GO TO 220
  200    IF (FF/FDEN.GE.9999.0D0) GO TO 240
  220    F(I) = FF/FDEN
         IFA = 1
         FP(I) = G01EDF('U',F(I),DBLE(IDF(I)),DBLE(IDF(3)),IFA)
  240 CONTINUE
      IFAIL = 0
      GO TO 280
  260 IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,0,P01REC)
  280 RETURN
      END
