      SUBROUTINE G13BDZ(R0,R,NL,NNA,S,NWDS,WA,NWA,WDS,ISF)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C        G13BDZ CARRIES OUT THE CALCULATIONS FOR G13BDF
C
C        R0     - CORRELATION AT LAG ZERO
C        R      - CORRELATION ARRAY
C        NL     - NO OF LAGS SUPPLIED
C        NNA    - ORDERS VECTOR
C        S      - STANDARD DEVIATION RATIO SY/SX
C        NWDS   - SIZE OF PARAMETER ARRAY
C        WA     - WORKING ARRAY
C        NWA    - SIZE OF WORKING ARRAY
C        WDS    - PARAMETER ARRAY
C        ISF    - SUCCESS/FAILURE PARAMETER
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  R0, S
      INTEGER           NL, NWA, NWDS
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NL), WA(NWA), WDS(NWDS)
      INTEGER           ISF(2), NNA(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  EF, EPS, PG, SM, X
      INTEGER           I, IB, IFAIL1, IP, IQ, J, K1, K2, K3, K4, K5,
     *                  KC, ME, ND
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F04ARF, G13AEX
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      DO 20 I = 1, NWDS
         WDS(I) = 10.0D0
   20 CONTINUE
C        SET UP SOME USEFUL VARIABLES
      IB = NNA(1)
      IQ = NNA(2)
      IP = NNA(3)
      ME = IP*IP
      ND = IQ + 2
      X = 1.0D0
      SM = X02AJF()
      EF = 1000.0D0
      EPS = EF*SM
C        INITIALISE ESTIMATION SUCCESS/FAILURE INDICATOR
      ISF(1) = 1
      ISF(2) = 0
      IF (IP.GT.0) ISF(2) = 1
C        ESTIMATE AR-LIKE PARAMETERS - THE DELTAS
      IF (IP.EQ.0) GO TO 200
C        EXPAND CORRELATIONS INTO TOEPLITZ MATRIX
      K1 = 0
      K2 = IB + IQ
      DO 100 I = 1, IP
         K3 = K2
         DO 80 J = 1, IP
            K1 = K1 + 1
            IF (K3.GE.IB) GO TO 40
            X = 0.0D0
            GO TO 60
   40       IF (K3.EQ.0) X = R0
            IF (K3.GT.0) X = R(K3)
   60       WA(K1) = X
            K3 = K3 + 1
   80    CONTINUE
         K2 = K2 - 1
  100 CONTINUE
C        MAKE UP RHS OF MATRIX EQUATION IN WDS(ND)
      K1 = ND
      K2 = IB + IQ
      DO 120 I = 1, IP
         K2 = K2 + 1
         WDS(K1) = R(K2)
         K1 = K1 + 1
  120 CONTINUE
C        SOLVE FOR DELTA PARAMETERS
      K1 = ME + 1
      IFAIL1 = 1
      CALL F04ARF(WA(1),IP,WDS(ND),IP,WDS(ND),WA(K1),IFAIL1)
      IF (IFAIL1.NE.0) GO TO 160
      K1 = ME + 1
      K2 = ND
C        CHECK DELTA PARAMETERS VALID
      DO 140 I = 1, IP
         WA(K1) = WDS(K2)
         K1 = K1 + 1
         K2 = K2 + 1
  140 CONTINUE
      K1 = ME + 1
      CALL G13AEX(WA(K1),IP,EPS,PG,KC)
      IF (KC.GE.0) GO TO 200
  160 ISF(2) = -1
      K1 = ND
      DO 180 I = 1, IP
         WDS(K1) = 0.0D0
         K1 = K1 + 1
  180 CONTINUE
  200 CONTINUE
C        ESTIMATE MA-LIKE PARAMETERS - THE OMEGAS
      IF (IB.EQ.0) X = R0
      IF (IB.GT.0) X = R(IB)
      WDS(1) = X*S
      IF (IQ.EQ.0) GO TO 280
      K1 = 1
      K2 = IB
      DO 240 I = 1, IQ
         K1 = K1 + 1
         K2 = K2 + 1
         WDS(K1) = R(K2)
         IF (ISF(2).LE.0) GO TO 240
         K5 = MIN(I,IP)
         K3 = K2
         K4 = ND
         DO 220 J = 1, K5
            K3 = K3 - 1
            IF (K3.EQ.0) X = R0
            IF (K3.GT.0) X = R(K3)
            WDS(K1) = WDS(K1) - WDS(K4)*X
            K4 = K4 + 1
  220    CONTINUE
  240 CONTINUE
      K1 = 1
      DO 260 I = 1, IQ
         K1 = K1 + 1
         WDS(K1) = -WDS(K1)*S
  260 CONTINUE
  280 CONTINUE
C        SUCCESSFUL EXIT
      RETURN
      END
