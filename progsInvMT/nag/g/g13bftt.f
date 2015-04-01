      SUBROUTINE G13BFT(MPAB,NMS,NUM,NUMW,MSN,MSPA,MSPB,NPAR,KEF,NAS,
     *                  NBS,NBVD,H,IH,WQ,IWQ,WZ,A,IDA,B,IDB,GC,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 16 REVISED. IER-1042 (JUN 1993).
C
C     SUBROUTINE G13BFT DERIVES THE CORRECTIONS TO THE
C     DERIVATIVE VECTOR ASSOCIATED WITH THE ARIMA PARAMETERS
C
C     .. Scalar Arguments ..
      INTEGER           IDA, IDB, IFAIL, IH, IWQ, KEF, NAS, NBS, NBVD,
     *                  NMS, NPAR, NUM, NUMW
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IDA), B(IDB), GC(NPAR), H(IH,NUM),
     *                  WQ(IWQ,NUM), WZ(NUM)
      INTEGER           MPAB(15), MSN(NMS), MSPA(NMS), MSPB(NMS)
C     .. Local Scalars ..
      DOUBLE PRECISION  DET, DETM, EPS, F, FC, GCQ, U, ZERO
      INTEGER           I, IRSPA, IRSPB, ISGN, ISN, IX, IXQ, J, JFA,
     *                  JFB, JSN, JX, JXQ, K, KFL, KFVM, KINC, KQ, KR,
     *                  KRSPA, KRSPB, KSGN, KSN, KT, KTFA, KTFB, KTSA,
     *                  KTSB, LAC, LACR, LBC, LBCR, LFA, LFB, LINC,
     *                  LRSPA, LRSPB, LSN, NUMR
C     .. Local Arrays ..
      INTEGER           MSNX(4)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F01ACZ, F03ABF, G13BFW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO/0.0D0/, U/1.0D0/
C     .. Executable Statements ..
C     KFVM IS THE SUBSCRIPT OF THE LAST NON-ARIMA PARAMETER
C     IN MSN,ETC
C
      KFVM = NMS - NPAR
      LAC = NAS - NBVD
      LBC = NBS - NBVD
      MSNX(1) = MPAB(2)
      MSNX(2) = -MPAB(3)
      MSNX(3) = MPAB(4)
      MSNX(4) = -MPAB(5)
      EPS = X02AJF()
C
C     CALCULATE DET, THE DETERMINANT OF H, FOR B.F.S ONLY
C     IF USING EXACT ESTIMATION AND FOR B.F.S+C+SIMPLE OMEGA
C     IF USING MARGINAL ESTIMATION
C
      CALL F03ABF(H,IH,NUM,DET,WZ,IFAIL)
      DETM = U/DET
C
C     THE ELEMENTS OF THE INVERSE OF H ARE HELD IN
C     ITS LOWER TRIANGLE
C
      CALL F01ACZ(NUM,EPS,H,IH,WQ,IWQ,WZ,I,IFAIL)
C
C     PROCESS EACH ARIMA PARAMETER IN TURN
C
      DO 420 K = 1, NPAR
         GC(K) = ZERO
         KQ = KFVM + K
C
C        DERIVE THE INCREMENT NEEDED TO OBTAIN THE A(T),B(T) SET
C        NUMBER FOR W*ARIMA COMBINATIONS
C
         DO 20 LINC = 1, 4
            IF (MSN(KQ).EQ.MSNX(LINC)) GO TO 40
   20    CONTINUE
         GO TO 80
   40    KINC = 0
         DO 60 I = 1, LINC
            IF (MSNX(I).EQ.0) GO TO 60
            KINC = KINC + 1
   60    CONTINUE
C
C        SUBSCRIPTS BEGINNING WITH K RELATE TO ARIMA PARAMETERS
C
   80    KRSPA = MSPA(KQ)
         KRSPB = MSPB(KQ)
         KSN = ABS(MSN(KQ))
         KSGN = 1
         IF (MSN(KQ).GE.0) GO TO 100
         KSGN = -1
C
C        PROCESS EACH ROW OF H IN TURN
C
  100    DO 400 I = 1, NUM
C
C           SUBSCRIPTS BEGINNING WITH I RELATE TO W,C, AND SIMPLE
C           OMEGA PARAMETERS
C
            IRSPA = MSPA(I+1)
            IRSPB = MSPB(I+1)
            ISN = ABS(MSN(I+1))
            ISGN = 1
            IF (MSN(I+1).GE.0) GO TO 120
            ISGN = -1
C
C           SUBSCRIPTS BEGINNING WITH L RELATE TO K*I COMBINATIONS
C
  120       LRSPA = KRSPA + IRSPA
            LRSPB = KRSPB + IRSPB
            LSN = ISN + KINC
            LFA = (LSN-1)*NAS + LRSPA + NBVD
            LFB = (LSN-1)*NBS + LRSPB + NBVD
C
C           PROCESS EACH COLUMN OF H IN TURN
C
            DO 380 J = 1, NUM
C
C              SUBSCRIPTS BEGINNING WITH J RELATE TO W,C AND SIMPLE
C              OMEGA PARAMETERS
C
               JSN = ABS(MSN(J+1))
               JFA = (JSN-1)*NAS + MSPA(J+1) + NBVD
               JFB = (JSN-1)*NBS + MSPB(J+1) + NBVD
               GCQ = ZERO
               NUMR = 1
C
C              PROCESS ALL J S IF I=1 AND ONLY J=1 OTHERWISE
C
               IF (I.GT.NUMW) GO TO 160
               IF (J.GT.NUMW) GO TO 160
               IF (I.NE.1) GO TO 140
               NUMR = NUMW + 1 - J
               GO TO 160
  140          IF (J.NE.1) GO TO 380
               NUMR = NUMW + 1 - I
C
C              CALCULATE F, THE RELEVANT ELEMENT OF THE INVERSE
C              OF MATRIX H
C
  160          CALL G13BFW(H,IH,I,J,KSGN,NUM,F,KFL,MSN,NMS,IX,JX)
               IF (F.NE.ZERO) GO TO 180
               IF (NUMR.EQ.1) GO TO 380
  180          LACR = LAC
               LBCR = LBC
               IF (NUMR.EQ.1) GO TO 200
               LACR = LAC - NUMR
               LBCR = LBC - NUMR
C
C              BUILD UP GCQ TO HOLD SUMS OF PRODUCTS OF A(T) COLUMNS
C              LESS SUMS OF PRODUCTS OF B(T) COLUMNS FOR THE
C              NUMBER OF VALUES IN THE A(T),B(T) SETS,LESS THE NUMBER
C              OF REPEATS
C
  200          DO 220 KT = 1, LACR
                  KTFA = JFA + KT
                  KTSA = LFA + KT
                  GCQ = GCQ + A(KTFA)*A(KTSA)
  220          CONTINUE
               KTFB = JFB + LBCR
               KTSB = LFB + LBCR
               IF (LBCR.LE.0) GO TO 260
               DO 240 KT = 1, LBCR
                  KTFB = JFB + KT
                  KTSB = LFB + KT
                  GCQ = GCQ - B(KTFB)*B(KTSB)
  240          CONTINUE
C
C              BUILD UP GCQ FURTHER TO TAKE ACCOUNT OF REPEATS
C              AND THEN DERIVE GC, THE OUTPUT CORRECTION VECTOR
C
  260          IF (NUMR.NE.1) GO TO 280
               GC(K) = GC(K) + F*GCQ
               GO TO 380
  280          FC = F
               DO 360 KR = 1, NUMR
                  KTFA = KTFA + 1
                  KTSA = KTSA + 1
                  GCQ = GCQ + A(KTFA)*A(KTSA)
                  KTFB = KTFB + 1
                  KTSB = KTSB + 1
                  LBCR = LBCR + 1
                  IF (LBCR.LE.0) GO TO 300
                  GCQ = GCQ - B(KTFB)*B(KTSB)
  300             IF (KR.NE.NUMR) GO TO 320
                  F = FC
                  GO TO 340
  320             IXQ = IX + NUMR - KR
                  JXQ = JX + NUMR - KR
                  F = H(IXQ,JXQ)
                  IF (KFL.EQ.1) GO TO 340
                  F = -F
  340             GC(K) = GC(K) + F*GCQ
  360          CONTINUE
  380       CONTINUE
  400    CONTINUE
  420 CONTINUE
      RETURN
      END
