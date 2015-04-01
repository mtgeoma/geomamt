      SUBROUTINE G13BFR(EX,ALPHA,A,PA,NA,IDA,W,BETA,B,IDW,IDB,WA,IDWA,
     *                  NRMP,KDQ,NP,NQ,NPS,NQS,NS,NPD,MPQS,NAS,NBS)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BFR CALCULATES THE DERIVATIVES W.R. TO
C     THE ARIMA PARAMETERS(I.E. THE A(T),B(T) SETS)
C
C     .. Scalar Arguments ..
      INTEGER           IDA, IDB, IDW, IDWA, KDQ, NA, NAS, NBS, NP, NPD,
     *                  NPS, NQ, NQS, NRMP, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IDA), ALPHA(NA), B(IDB), BETA(IDW), EX(NA),
     *                  PA(NA), W(IDW), WA(IDWA)
      INTEGER           MPQS(4)
C     .. Local Scalars ..
      INTEGER           I, ID, JA, JB, KCT, KWSPH, KWSTH, KWTH
C     .. External Subroutines ..
      EXTERNAL          G13AEU
C     .. Executable Statements ..
      KCT = 0
      KWTH = 1 + NRMP
      KWSPH = KWTH + NRMP
      KWSTH = KWSPH + NRMP
C
C     PROCESS EACH OF THE FOUR TYPES OF ARIMA PARAMETER
C     IN TURN
C
      DO 40 I = 1, 4
         IF (MPQS(I).LE.0) GO TO 40
         KCT = KCT + 1
         JA = 1 + (KCT-1)*NAS + KDQ
         JB = 1 + (KCT-1)*NBS + KDQ
         ID = 2 + I
C
C        CALCULATE THE A(T),B(T) SETS IN ONE WAY FOR AUTOREGRESSIVE
C        PARAMETERS AND IN ANOTHER WAY FOR MOVING AVERAGE PARAMETERS
C
         IF (ID.EQ.(2*(ID/2))) GO TO 20
         CALL G13AEU(ID,EX,ALPHA,A(JA),NA,W,BETA,B(JB),IDW,WA(1)
     *               ,WA(KWTH),WA(KWSPH),WA(KWSTH)
     *               ,NRMP,NP,NQ,NPS,NQS,NS,NPD)
         GO TO 40
   20    CALL G13AEU(ID,PA,ALPHA,A(JA),NA,W,BETA,B(JB),IDW,WA(1)
     *               ,WA(KWTH),WA(KWSPH),WA(KWSTH)
     *               ,NRMP,NP,NQ,NPS,NQS,NS,NPD)
   40 CONTINUE
      RETURN
      END
