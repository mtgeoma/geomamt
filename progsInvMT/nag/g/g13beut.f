      SUBROUTINE G13BEU(LAC,LBC,MSN,MASPA,MASPB,NMSU,A,IDA,B,IDB,NBVD,S,
     *                  G,IGH,H,IH,MPAB)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEU CALCULATES S, G, AND H USING THE
C     RELEVANT A(T) AND B(T) SETS, AND THEIR START POINTS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S
      INTEGER           IDA, IDB, IGH, IH, LAC, LBC, NBVD, NMSU
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IDA), B(IDB), G(IGH), H(IH,IGH)
      INTEGER           MASPA(IH), MASPB(IH), MPAB(15), MSN(IH)
C     .. Local Scalars ..
      DOUBLE PRECISION  T, U, ZERO
      INTEGER           I, IPHA, IPHB, IQ, ISPHA, ISPHB, J, JSA, JSB,
     *                  JSQA, JSQB, K, KK, KPH, KSA, KSB, KSPH, KSQA,
     *                  KSQB, NAS, NBS
C     .. Data statements ..
      DATA              ZERO/0.0D0/, U/1.0D0/
C     .. Executable Statements ..
C
C     NAS AND NBS CONTAIN THE FULL LENGTHS OF EACH A(T)
C     AND B(T) SET HELD IN A AND B.
C
      NAS = LAC + NBVD
      NBS = LBC + NBVD
C
C     ZEROISE S, G, AND H
C
      S = ZERO
      KK = NMSU - 1
      IF (NMSU.LE.1) GO TO 60
      DO 40 I = 1, KK
         G(I) = ZERO
         DO 20 J = 1, KK
            H(J,I) = ZERO
   20    CONTINUE
   40 CONTINUE
C
C     FORM S
C
   60 DO 80 I = 1, LAC
         IQ = MASPA(1) - 1 + I
         S = S + A(IQ)*A(IQ)
   80 CONTINUE
      IF (LBC.LE.0) GO TO 120
      DO 100 I = 1, LBC
         IQ = MASPB(1) - 1 + I
         S = S - B(IQ)*B(IQ)
  100 CONTINUE
  120 IF (NMSU.LE.1) GO TO 380
C
C     PROCESS NMSU SETS OF A(T) AND B(T) AND START POINT COMBINATIONS
C
      DO 340 J = 1, NMSU
         JSA = MASPA(J)
         JSB = MASPB(J)
C
C        AFTER THE S CALCULATION IT SETS THE FIRST NBS VALUES
C        OF A(T) AND B(T) TO ZERO WHERE PHI AND SPHI ARE INVOLVED
C
         IF (J.NE.2) GO TO 180
         IF (NBS.LE.0) GO TO 180
         KPH = MPAB(2)
         IPHA = (KPH-1)*NAS
         IPHB = (KPH-1)*NBS
         KSPH = MPAB(4)
         ISPHA = (KSPH-1)*NAS
         ISPHB = (KSPH-1)*NBS
         IF ((KPH+KSPH).LE.0) GO TO 180
         DO 160 I = 1, NBS
            IF (KPH.LE.0) GO TO 140
            IPHA = IPHA + 1
            IPHB = IPHB + 1
            A(IPHA) = ZERO
            B(IPHB) = ZERO
  140       IF (KSPH.LE.0) GO TO 160
            ISPHA = ISPHA + 1
            ISPHB = ISPHB + 1
            A(ISPHA) = ZERO
            B(ISPHB) = ZERO
  160    CONTINUE
  180    DO 320 K = J, NMSU
            IF (K.EQ.1) GO TO 320
            KSA = MASPA(K)
            KSB = MASPB(K)
C
C           FORM SUMS OF SQUARES OR PRODUCTS
C
            T = ZERO
            DO 200 I = 1, LAC
               JSQA = JSA - 1 + I
               KSQA = KSA - 1 + I
               T = T + A(JSQA)*A(KSQA)
  200       CONTINUE
            IF (LBC.LE.0) GO TO 240
            DO 220 I = 1, LBC
               JSQB = JSB - 1 + I
               KSQB = KSB - 1 + I
               T = T - B(JSQB)*B(KSQB)
  220       CONTINUE
C
C           CHANGE SIGN OF T IF ONE OF THE MSN VALUES IS ZERO
C
  240       IF (MSN(J).GT.0) GO TO 260
            T = -T
  260       IF (MSN(K).GT.0) GO TO 280
            T = -T
  280       IF (J.NE.1) GO TO 300
C
C           G HOLDS SP OF Y WITH EACH X
C
            G(K-1) = T
            GO TO 320
C
C           H HOLDS SP OF EACH PAIR OF XS
C
  300       H(J-1,K-1) = T
  320    CONTINUE
  340 CONTINUE
      DO 360 I = 1, KK
         IF (H(I,I).NE.ZERO) GO TO 360
         H(I,I) = U
  360 CONTINUE
  380 RETURN
      END
