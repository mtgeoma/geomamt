      SUBROUTINE G13AEW(EX,ALPHA,NA,BETA,NB,AQ,NAQ,BQ,NBQ,KSCH,PHI,
     *                  THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD,NQD)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AEW CALCULATES THE A S AND B S NEEDED TO CALCULATE THE
C     DERIVATIVES OF S. THESE ARE HELD IN THE AQ AND BQ ARRAYS
C     AND ARE SUBDIVIDED INTO SEVEN SECTIONS RELATING TO THE X S,
C     THE PRE-X S, THE PHI S, THE THETA S, THE SPHI S, THE STHETA S
C     AND THE CONSTANT. THOSE SECTIONS NOT NEEDED ARE STORED
C     AS ZERO S
C
C     .. Scalar Arguments ..
      INTEGER           KSCH, NA, NAQ, NB, NBQ, NP, NPD, NPS, NQ, NQD,
     *                  NQS, NRMP, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  ALPHA(NA), AQ(NAQ), BETA(NB), BQ(NBQ), EX(NA),
     *                  PHI(NRMP), SPHI(NRMP), STHETA(NRMP), THETA(NRMP)
C     .. Local Scalars ..
      DOUBLE PRECISION  U, ZERO
      INTEGER           I, J, K, NX
C     .. External Subroutines ..
      EXTERNAL          G13AEU
C     .. Data statements ..
      DATA              ZERO/0.0D0/, U/1.0D0/
C     .. Executable Statements ..
C
C     ZEROISE AQ AND BQ
C
      DO 20 I = 1, NAQ
         AQ(I) = ZERO
   20 CONTINUE
      DO 40 I = 1, NBQ
         BQ(I) = ZERO
   40 CONTINUE
C
C     USE THE EXTENDED SERIES EX TO OBTAIN A AND B
C
      CALL G13AEU(1,EX,ALPHA,AQ(1),NA,BQ(1),BETA,BQ(1)
     *            ,NB,PHI,THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
      IF (KSCH.EQ.1) GO TO 200
      IF (NQD.LE.0) GO TO 60
C
C     USE A SPECIAL SERIES TO OBTAIN AX AND BX
C
      AQ(2*NA+1) = U
      CALL G13AEU(2,AQ(2*NA+1),ALPHA,AQ(NA+1),NA,BQ(1),BETA,BQ(NB+1)
     *            ,NB,PHI,THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
      AQ(2*NA+1) = ZERO
   60 IF (KSCH.EQ.2) GO TO 200
      IF (KSCH.EQ.3) GO TO 120
      K = NQD + 2*NA
      NX = NA - NQD
      DO 80 I = 1, NX
         J = I + K
         AQ(J) = U
   80 CONTINUE
C
C     USE A SPECIAL SERIES TO OBTAIN AC AND BC
C
      CALL G13AEU(7,AQ(2*NA+1),ALPHA,AQ(6*NA+1),NA,BQ(1),BETA,BQ(6*NB+1)
     *            ,NB,PHI,THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
      K = 2*NA
      DO 100 I = 1, NA
         J = I + K
         AQ(J) = ZERO
  100 CONTINUE
  120 IF (NP.LE.0) GO TO 140
C
C     USE THE EXTENDED SERIES EX TO OBTAIN THE A S AND B S RELATING
C     TO THE PHI S
C
      CALL G13AEU(3,EX,ALPHA,AQ(2*NA+1),NA,BQ(1),BETA,BQ(2*NB+1)
     *            ,NB,PHI,THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
  140 IF (NQ.LE.0) GO TO 160
C
C     USE A AND B TO OBTAIN THE A S AND B S RELATING TO THE THETA S
C
      CALL G13AEU(4,AQ(1),ALPHA,AQ(3*NA+1),NA,BQ(1),BETA,BQ(3*NB+1)
     *            ,NB,PHI,THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
  160 IF (NPS.LE.0) GO TO 180
C
C     USE THE EXTENDED SERIES EX TO OBTAIN THE A S AND B S
C     RELATING TO THE SPHI S
C
      CALL G13AEU(5,EX,ALPHA,AQ(4*NA+1),NA,BQ(1),BETA,BQ(4*NB+1)
     *            ,NB,PHI,THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
  180 IF (NQS.LE.0) GO TO 200
C
C     USE A AND B TO OBTAIN THE A S AND B S RELATING TO THE STHETA S
C
      CALL G13AEU(6,AQ(1),ALPHA,AQ(5*NA+1),NA,BQ(1),BETA,BQ(5*NB+1)
     *            ,NB,PHI,THETA,SPHI,STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
  200 RETURN
      END
