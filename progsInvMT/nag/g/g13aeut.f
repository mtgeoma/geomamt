      SUBROUTINE G13AEU(ID,EX,ALPHA,A,NA,W,BETA,B,NB,PHI,THETA,SPHI,
     *                  STHETA,NRMP,NP,NQ,NPS,NQS,NS,NPD)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11 REVISED. IER-440 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-520 (AUG 1986).
C
C     G13AEU CALCULATES A AND B ARRAYS RELATING TO THE DIFFERENT
C     TYPES OF PARAMETER SETS
C
C     .. Scalar Arguments ..
      INTEGER           ID, NA, NB, NP, NPD, NPS, NQ, NQS, NRMP, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NA), ALPHA(NA), B(*), BETA(NB), EX(NA),
     *                  PHI(NRMP), SPHI(NRMP), STHETA(NRMP),
     *                  THETA(NRMP), W(NB)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
      INTEGER           I, J, KA, KALPHA, KB, KBETA, KQ, KX
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
C
C     CALCULATE ALPHA(I) VALUES, FOR I=1 TO NA, WHERE NA IS THE
C     NUMBER OF PRE-X S PLUS THE NUMBER OF X S
C
      DO 80 I = 1, NA
C
C        OBTAIN INITIAL VALUE OF ALPHA(I)
C
         ALPHA(I) = EX(I)
         IF (ID.EQ.5) GO TO 40
         IF (ID.EQ.6) GO TO 40
         IF (ID.EQ.4) GO TO 80
         IF (NPS.LE.0) GO TO 40
C
C        BUILD UP ALPHA(I) USING SPHI PARAMETER VALUES
C
         DO 20 J = 1, NPS
            KX = I - J*NS
            IF (KX.LE.0) GO TO 40
            ALPHA(I) = ALPHA(I) - SPHI(J)*EX(KX)
   20    CONTINUE
   40    IF (NQS.LE.0) GO TO 80
C
C        BUILD UP ALPHA(I) USING STHETA PARAMETER VALUES
C
         DO 60 J = 1, NQS
            KALPHA = I - J*NS
            IF (KALPHA.LE.0) GO TO 80
            ALPHA(I) = ALPHA(I) + STHETA(J)*ALPHA(KALPHA)
   60    CONTINUE
   80 CONTINUE
C
C     CALCULATE A(I) VALUES, FOR I=1 TO NA
C
      DO 160 I = 1, NA
C
C        OBTAIN INITIAL VALUE OF A(I)
C
         A(I) = ALPHA(I)
         IF (ID.EQ.4) GO TO 120
         IF (ID.EQ.3) GO TO 120
         IF (ID.EQ.6) GO TO 160
         IF (NP.LE.0) GO TO 120
C
C        BUILD UP A(I) USING PHI PARAMETER VALUES
C
         DO 100 J = 1, NP
            KALPHA = I - J
            IF (KALPHA.LE.0) GO TO 120
            A(I) = A(I) - PHI(J)*ALPHA(KALPHA)
  100    CONTINUE
  120    IF (NQ.LE.0) GO TO 160
C
C        BUILD UP A(I) USING THETA PARAMETER VALUES
C
         DO 140 J = 1, NQ
            KA = I - J
            IF (KA.LE.0) GO TO 160
            A(I) = A(I) + THETA(J)*A(KA)
  140    CONTINUE
  160 CONTINUE
      IF (ID.LE.0) GO TO 600
      IF (NPD.LE.0) GO TO 600
      KQ = NPS*NS
C
C     CALCULATE BETA(I) VALUES, FOR I=1 TO NPD, WHERE NPD IS THE
C     LENGTH OF THE BETA AND B SETS
C
      DO 380 I = 1, NPD
C
C        OBTAIN INITIAL VALUE OF BETA(I)
C
         IF (ID.EQ.7) GO TO 240
         IF (ID-4) 240, 180, 180
  180    IF (ID-5) 200, 220, 200
  200    BETA(I) = W(I)
         IF (ID-6) 380, 340, 380
  220    BETA(I) = EX(I)
         GO TO 340
  240    KX = I - KQ
         IF (KX) 260, 260, 280
  260    BETA(I) = ZERO
         GO TO 300
  280    BETA(I) = EX(KX)
  300    IF (NPS.LE.0) GO TO 340
C
C        BUILD UP BETA(I) USING SPHI PARAMETER VALUES
C
         DO 320 J = 1, NPS
            KX = I + J*NS - KQ
            IF (KX.LE.0) GO TO 320
            BETA(I) = BETA(I) - SPHI(J)*EX(KX)
  320    CONTINUE
  340    IF (NQS.LE.0) GO TO 380
C
C        BUILD UP BETA(I) USING STHETA PARAMETER VALUES
C
         DO 360 J = 1, NQS
            KBETA = I - J*NS
            IF (KBETA.LE.0) GO TO 380
            BETA(I) = BETA(I) + STHETA(J)*BETA(KBETA)
  360    CONTINUE
  380 CONTINUE
C
C     CALCULATE B(I)VALUES, FOR I=1 TO NPD
C
      DO 580 I = 1, NPD
C
C        OBTAIN INITIAL VALUE OF B(I)
C
         IF (ID.EQ.7) GO TO 440
         IF (ID-3) 440, 400, 400
  400    IF (ID-5) 420, 440, 420
  420    B(I) = BETA(I)
         IF (ID-6) 540, 580, 540
  440    KBETA = I - NP
         IF (KBETA) 460, 460, 480
  460    B(I) = ZERO
         GO TO 500
  480    B(I) = BETA(KBETA)
  500    IF (NP.LE.0) GO TO 540
C
C        BUILD UP B(I) USING PHI PARAMETER VALUES
C
         DO 520 J = 1, NP
            KBETA = I + J - NP
            IF (KBETA.LE.0) GO TO 520
            B(I) = B(I) - PHI(J)*BETA(KBETA)
  520    CONTINUE
  540    IF (NQ.LE.0) GO TO 580
C
C        BUILD UP B(I) USING THETA PARAMETER VALUES
C
         DO 560 J = 1, NQ
            KB = I - J
            IF (KB.LE.0) GO TO 580
            B(I) = B(I) + THETA(J)*B(KB)
  560    CONTINUE
  580 CONTINUE
  600 RETURN
      END
