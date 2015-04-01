      SUBROUTINE G13ADZ(R,NL,MR,YV,PAR,NPAR,WA,NWA,RV,ISF)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     CALCULATIONS FOR G13ADF
C
C     PARAMETERS
C     R       - ARRAY OF AUTOCORRELATIONS
C     NL      - NO. OF AUTOCORRELATIONS
C     MR      - ORDERS VECTOR, DIMENSION 7
C     YV      - SAMPLE VARIANCE
C     PAR     - PARAMETERS VECTOR
C     NPAR    - EXACT NUMBER OF PARAMETERS
C     WA      - WORKING ARRAY
C     NWA     - SIZE OF WORKING ARRAY
C     RV      - RESIDUAL VARIANCE
C     ISF     - ESTIMATION SUCCESS INDICATOR ARRAY
C
C     USES NAG LIBRARY ROUTINES F03AFF, F04AJF, G13ADY, G13AEX,
C     X02AJF
C
C
C     INITIALISE ESTIMATION SUCCESS INDICATORS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RV, YV
      INTEGER           NL, NPAR, NWA
C     .. Array Arguments ..
      DOUBLE PRECISION  PAR(NPAR), R(NL), WA(NWA)
      INTEGER           ISF(4), MR(7)
C     .. Local Scalars ..
      DOUBLE PRECISION  DET, EF, EPS, EPSILN, PG, Q, SM
      INTEGER           I, I1, I2, ID, IFAIL1, INDK, INDPAR, INDR,
     *                  INDTP, ISEAS, ISP, ISQ, IWA, J1, J2, J3, JWA,
     *                  KC, KWA, L1, L2, MAXITN
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F03AFF, F04AJF, G13ADY, G13AEX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      DO 20 I = 1, 4
         ISF(I) = 0
   20 CONTINUE
      IF (MR(1).GT.0) ISF(1) = 1
      IF (MR(3).GT.0) ISF(2) = 1
      IF (MR(4).GT.0) ISF(3) = 1
      IF (MR(6).GT.0) ISF(4) = 1
C
C     FURTHER INITIALISATION
C
      L1 = 0
      L2 = 0
      ISEAS = 1
      Q = 1.0D0
      SM = X02AJF()
      EF = 1000.0D0
      EPS = EF*SM
      RV = YV
C
C     LOOP TWICE - NON-SEASONAL,SEASONAL ESTIMATION
C
      DO 600 I = 1, 2
C
C        CALC MR INDICES AND SEASONALITY
C
         ISP = MR(3*I-2)
         ISQ = MR(3*I)
         IF (I.EQ.2) ISEAS = MR(7)
C
C        JUMP IF NO AUTOREGRESSIVE PARAMETERS
C
         IF (ISP.EQ.0) GO TO 240
C
C        UPDATE PAR ARRAY POINTERS
C
         L1 = L2 + 1
         L2 = L1 + ISP - 1
C
C        ESTIMATE AUTOREGRESSIVE PARAMETERS
C
         DO 40 I1 = L1, L2
            PAR(I1) = 0.0D0
   40    CONTINUE
C
C        EXPAND CORRELATIONS INTO TOEPLITZ MATRIX
C
         INDTP = 0
         DO 140 I1 = 1, ISP
            DO 120 I2 = 1, ISP
               INDTP = INDTP + 1
               INDR = (ISQ-I1+I2)*ISEAS
               IF (INDR) 80, 60, 80
   60          WA(INDTP) = 1.0D0
               GO TO 100
   80          INDR = ABS(INDR)
               WA(INDTP) = R(INDR)
  100          CONTINUE
  120       CONTINUE
  140    CONTINUE
C
C        EXPAND CORRS INTO PAR ARRAY AT L1
C
         INDPAR = L1
         DO 160 I1 = 1, ISP
            INDR = (ISQ+I1)*ISEAS
            PAR(INDPAR) = R(INDR)
            INDPAR = INDPAR + 1
  160    CONTINUE
C
C        SOLVE FOR AUTOREGRESSIVE PARAMETERS
C
         IWA = INDTP + 1
         JWA = INDTP + ISP
         IFAIL1 = 1
         CALL F03AFF(ISP,SM,WA(1),ISP,DET,ID,WA(IWA),IFAIL1)
         IF (IFAIL1.NE.0) GO TO 200
         CALL F04AJF(ISP,1,WA(1),ISP,WA(IWA),PAR(L1),ISP)
C        CHECK AUTOREGRESSIVE PARAMETERS VALID
         INDPAR = L1
         DO 180 I1 = IWA, JWA
            WA(I1) = PAR(INDPAR)
            INDPAR = INDPAR + 1
  180    CONTINUE
         CALL G13AEX(WA(IWA),ISP,EPS,PG,KC)
         IF (KC.GE.0) GO TO 240
  200    ISF(2*I-1) = -1
         DO 220 I1 = L1, L2
            PAR(I1) = 0.0D0
  220    CONTINUE
  240    CONTINUE
C
C        MODIFY AUTOCORRELATIONS FOR AUTOREGRESSIVE PARAMETERS
         INDTP = ISQ + ISP + 1
         DO 440 J1 = 1, INDTP
            INDK = J1 - 1
            IF (INDK-ISQ) 280, 280, 260
  260       WA(J1) = 0.0D0
            GO TO 440
  280       IF (INDK) 300, 300, 320
  300       WA(J1) = 1.0D0
            GO TO 340
  320       INDK = INDK*ISEAS
            WA(J1) = R(INDK)
  340       CONTINUE
            IF (ISF(2*I-1).LE.0) GO TO 440
            INDPAR = L1
            DO 420 J2 = 1, ISP
               INDR = INDK - (J2*ISEAS)
               IF (INDR) 380, 360, 380
  360          WA(J1) = WA(J1) - PAR(INDPAR)
               GO TO 400
  380          INDR = ABS(INDR)
               WA(J1) = WA(J1) - PAR(INDPAR)*R(INDR)
  400          INDPAR = INDPAR + 1
  420       CONTINUE
  440    CONTINUE
         INDTP = ISQ + 1
         IF (ISF(2*I-1).LE.0) GO TO 500
         DO 480 J1 = 1, INDTP
            INDPAR = L1
            DO 460 J2 = 1, ISP
               J3 = J1 + J2
               IF ((J3).GT.INDTP) GO TO 480
               WA(J1) = WA(J1) - PAR(INDPAR)*WA(J3)
               INDPAR = INDPAR + 1
  460       CONTINUE
  480    CONTINUE
  500    CONTINUE
C
C        MODIFY RV IF THERE ARE NO M.A. PARAMETERS
C
         IF (ISQ.EQ.0) RV = RV*WA(1)
C
C        JUMP IF NO MOVING AVERAGE PARAMETERS
C
         IF (ISQ.EQ.0) GO TO 580
C
C        UPDATE PAR ARRAY POINTERS
C
         L1 = L2 + 1
         L2 = L1 + ISQ - 1
C
C        ESTIMATE MOVING AVERAGE PARAMETERS
C
C        G13ADY RETURNS THE UNSCALED RESIDUAL VARIANCE AND MOVING
C        AVERAGE PARAMETERS (RESP.) IN WA(1) TO WA(ISQ+1)
C
         IWA = INDTP + 1
         JWA = IWA + INDTP
         KWA = JWA + INDTP
         EPSILN = MAX(100.0D0*X02AJF(),1.0D-7)
         MAXITN = 20
         CALL G13ADY(WA(1),WA(IWA),WA(JWA),WA(KWA)
     *               ,INDTP,EPSILN,MAXITN,IFAIL1)
         RV = RV*WA(1)
         IF (IFAIL1.NE.0) GO TO 540
C
C        TRANSFER PARAMETERS TO PAR
C
         INDPAR = L1
         DO 520 J1 = 1, ISQ
            PAR(INDPAR) = WA(J1+1)
            INDPAR = INDPAR + 1
  520    CONTINUE
C
C        CHECK MOVING AVERAGE PARAMETERS VALID
C
         CALL G13AEX(WA(2),ISQ,EPS,PG,KC)
         IF (KC.GT.0) GO TO 580
  540    ISF(2*I) = -1
         INDPAR = L1
         DO 560 J1 = 1, ISQ
            PAR(INDPAR) = 0.0D0
            INDPAR = INDPAR + 1
  560    CONTINUE
  580    CONTINUE
  600 CONTINUE
      RETURN
      END
