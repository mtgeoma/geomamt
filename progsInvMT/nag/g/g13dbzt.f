      SUBROUTINE G13DBZ(C0,LC0,C,LC,NSM,NS,NL,NK,P,V0,V,D,LW,DB,W,WB,
     *                  NVP,WA,NWA,IERROR)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        G13DBZ CARRIES OUT THE CALCULATIONS FOR G13DBF
C
C        C0     - LAG ZERO COVARIANCES
C        C      - CROSS COVARIANCE MATRIX, LAGS 1 TO NL
C        NSM    - MAX. NO. OF SERIES
C        NS     - NO. OF SERIES
C        NL     - NO. OF LAGS FOR CROSS COVS
C        NK     - NO. OF LAGS FOR PARTIAL CORRS.
C        P      - MULTIVARIATE PARTIAL AUTOCORRS
C        V0     - LAG 0 PRED. ERR. VARIANCE DET.
C        V      - PRED. ERR. VARIANCE RATIOS
C        D      - PRED. ERR. VARIANCES
C        DB     - LAST BACK PRED. ERR. VARIANCES
C        W      - PREDICTION COEFFICIENTS
C        WB     - BACKWARDS PREDICTION COEFFICIENTS
C        NVP    - NO. OF VALID PARAMETERS
C        WA     - WORK ARRAY
C        NWA    - SIZE OF WORK ARRAY
C        IERROR - ERROR INDICATOR
C
C        INITIALIZE LOCAL VARIABLES
C     .. Scalar Arguments ..
      DOUBLE PRECISION  V0
      INTEGER           IERROR, LC, LC0, LW, NK, NL, NS, NSM, NVP, NWA
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LC), C0(LC0), D(LW), DB(LC0), P(NK), V(NK),
     *                  W(LW), WA(NWA), WB(LW)
C     .. Local Scalars ..
      DOUBLE PRECISION  DET, DETL
      INTEGER           IFAIL1, J, K, K1, K2, K3, NK1, NSMQ, NSMQK, NSQ,
     *                  NSQ1, NWA1
C     .. External Subroutines ..
      EXTERNAL          F03ABF, G13DBU, G13DBV, G13DBW, G13DBX, G13DBY
C     .. Executable Statements ..
      NSMQ = NSM*NSM
      NSQ = NS*NS
      NSQ1 = NSQ + 1
      NWA1 = NWA - NSQ
      NK1 = NK - 1
      NVP = 0
C        ** STEP 0 OF RECURSION **
C        COPY DBAR (=C ) INTO DB
C                 0   0
      CALL G13DBV(C0,NSM,DB,NSM,NS)
C        CALCULATE V0
      IFAIL1 = 1
      CALL F03ABF(DB,NSM,NS,DET,WA,IFAIL1)
C        RESTORE DB TO SYMMETRY
      CALL G13DBV(DB,NSM,DB,NSM,NS)
      IF (IFAIL1.EQ.0) GO TO 40
      IERROR = 2
C        SET MATRICES AND VECTORS TO ZERO
      K = 1
      CALL G13DBU(P,V,D,W,WB,K,NSM,NS,NK)
C        SET V0 AND DB TO ZERO
      V0 = 0.0D0
      DO 20 J = 1, NSQ
         DB(J) = 0.0D0
   20 CONTINUE
      GO TO 180
   40 V0 = DET
C        ** STEP 1 OF RECURSION **
      K = 1
C        CREATE G (=C ) IN WORK STORE
C                1   1
      CALL G13DBY(C,W,WA,NSM,NS,K,WA(NSQ1),NWA1,NSMQ,NSQ)
C        CALCULATE PHI, PHIBAR, D, AND DBAR
      K2 = 1
      K3 = K2 + NSMQ
      CALL G13DBW(W(K2),WB(K2),DB(1),DB(1),D(K2),D(K3),NSM,WA(1)
     *            ,NSM,NS,WA(NSQ1),NWA1)
C        CHECK D  POSITIVE DEFINITE
C               1
      IFAIL1 = 1
      CALL F03ABF(D(K2),NSM,NS,DET,WA,IFAIL1)
      IF (IFAIL1.EQ.0) GO TO 60
C        IF NOT, SET MATRICES AND VECTORS TO ZERO
      CALL G13DBU(P,V,D,W,WB,K,NSM,NS,NK)
      GO TO 160
   60 CONTINUE
C        IF SO, COPY DBAR  INTO DB
C                        1
      CALL G13DBV(D(K3),NSM,DB,NSM,NS)
C        CALCULATE V AND P
      V(1) = DET/V0
      P(1) = 1.0D0 - DET/V0
      DETL = DET
C        INCREMENT NVP
      NVP = NVP + 1
C        ** STEPS 2 TO NK-1 OF RECURSION **
      IF (NK.LT.3) GO TO 120
      DO 100 K = 2, NK1
C        CREATE G  IN WORK SPACE
C                K
         NSMQK = NSMQ*K
         CALL G13DBY(C,W,WA,NSM,NS,K,WA(NSQ1),NWA1,NSMQK,NSQ)
C        CALCULATE PHI, PHIBAR, D, AND DBAR
         K1 = K2
         K2 = K3
         K3 = K3 + NSMQ
         CALL G13DBW(W(K2),WB(K2),D(K1),D(K2),D(K2),D(K3),NSM,WA(1)
     *               ,NSM,NS,WA(NSQ1),NWA1)
C        CHECK D  POSITIVE DEFINITE
C               K
         IFAIL1 = 1
         CALL F03ABF(D(K2),NSM,NS,DET,WA,IFAIL1)
         IF (IFAIL1.EQ.0) GO TO 80
C        IF NOT, SET MATRICES AND VECTORS TO ZERO
         CALL G13DBU(P,V,D,W,WB,K,NSM,NS,NK)
         GO TO 160
   80    CONTINUE
C        IF SO, COPY DBAR  INTO DB
C                        K
         CALL G13DBV(D(K3),NSM,DB,NSM,NS)
C        UPDATE PHI AND PHIBAR
         CALL G13DBX(W,WB,NSM,NS,K,WA,NWA,NSMQK)
C        CALCULATE V AND P
         V(K) = DET/V0
         P(K) = 1.0D0 - DET/DETL
         DETL = DET
C        INCREMENT NVP
         NVP = NVP + 1
  100 CONTINUE
  120 CONTINUE
C        ** STEP NK OF RECURSION **
      K = NK
      IF (NK.LT.2) GO TO 160
C        CREATE G   IN WORK SPACE
C                NK
      NSMQK = NSMQ*K
      CALL G13DBY(C,W,WA,NSM,NS,NK,WA(NSQ1),NWA1,NSMQK,NSQ)
C        CALCULATE PHI, PHIBAR, D, AND DBAR
      K1 = K2
      K2 = K3
      CALL G13DBW(W(K2),WB(K2),D(K1),D(K2),D(K2),WA(NSQ1),NS,WA(1)
     *            ,NSM,NS,WA(NSQ1),NWA1)
C        CHECK D   POSITIVE DEFINITE
C               NK
      IFAIL1 = 1
      CALL F03ABF(D(K2),NSM,NS,DET,WA,IFAIL1)
      IF (IFAIL1.EQ.0) GO TO 140
C        IF NOT, SET MATRICES AND VECTORS TO ZERO
      CALL G13DBU(P,V,D,W,WB,K,NSM,NS,NK)
      GO TO 160
  140 CONTINUE
C        IF SO, COPY DBAR   INTO DB
C                        NK
      CALL G13DBV(WA(NSQ1),NS,DB,NSM,NS)
C        UPDATE PHI AND PHIBAR
      NSMQK = NSMQ*K
      CALL G13DBX(W,WB,NSM,NS,K,WA,NWA,NSMQK)
C        RESTORE D   TO SYMMETRY
C                 NK
      CALL G13DBV(D(K2),NSM,D(K2),NSM,NS)
C        CALCULATE V AND P
      V(NK) = DET/V0
      P(NK) = 1.0D0 - DET/DETL
C        INCREMENT NVP
      NVP = NVP + 1
  160 CONTINUE
C        ERROR FOR PREMATURE TERMINATION
      IF (NVP.LT.NK) IERROR = 3
  180 CONTINUE
      RETURN
      END
