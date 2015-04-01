      SUBROUTINE G02HAZ(INDW,IPSI,INDC,SIGMA,N,M,X,IX,RS,WGT,CPSI,H1,H2,
     *                  H3,C,IC,D,E,WK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1033 (JUN 1993).
C     MARK 17 REVISED. IER-1677 (JUN 1995).
C
C     Calculate asymptotic variance-covariance matrix.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CPSI, H1, H2, H3, SIGMA
      INTEGER           IC, IFAIL, INDC, INDW, IPSI, IX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(IC,M), D(N), E(N), RS(N), WGT(N),
     *                  WK(M*(N+M+1)), X(IX,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, FACT, PS, RN, S, S1, S2, SUM2, TMP1, TMP2,
     *                  VAR, XKAPPA, XMU, XMU2, Z
      INTEGER           I, IFAIL2, IJ, J, J1, M1, MM1, NUMC
C     .. External Functions ..
      DOUBLE PRECISION  G02HAU, G07DBX, DDOT, X02AJF
      EXTERNAL          G02HAU, G07DBX, DDOT, X02AJF
C     .. External Subroutines ..
      EXTERNAL          F01AAZ, F01ACZ, F06FCF, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     SET UP PARAMETER VALUES
C
      IFAIL = 0
      IFAIL2 = 1
      EPS = X02AJF()
      M1 = M + 1
      MM1 = M1*M + 1
C
C     FOR HUBER TYPE REGRESSION
C
      IF (INDW.EQ.0) THEN
C
C        CALCULATE INV(X'X)
C
         DO 40 J = 1, M
            DO 20 I = 1, J
               IJ = (J-1)*(M+1) + I
               WK(IJ) = DDOT(N,X(1,I),1,X(1,J),1)
   20       CONTINUE
   40    CONTINUE
         CALL F01ACZ(M,EPS,WK,M1,WK(MM1),N,D,NUMC,IFAIL2)
         IF (IFAIL2.NE.0) THEN
            IFAIL = 1
            RETURN
         END IF
C
C        -------
C        COMPUTES CORRECTION FACTORS XKAPPA AND SUM2 FOR
C        THE COVARIANCE MATRIX.
C        FACS CALLS THE FUNCTIONS PSI AND PSP
C
         TMP1 = 0.0D0
         TMP2 = 0.0D0
         DO 60 J = 1, N
            S = RS(J)/SIGMA
            TMP1 = G02HAU(S,IPSI,CPSI,H1,H2,H3) + TMP1
            PS = G07DBX(S,IPSI,CPSI,H1,H2,H3)
            TMP2 = PS*PS + TMP2
   60    CONTINUE
         XMU = TMP1/DBLE(N)
         SUM2 = TMP2
         VAR = 0.0D0
         DO 80 J = 1, N
            S = RS(J)/SIGMA
            VAR = (G02HAU(S,IPSI,CPSI,H1,H2,H3)-XMU)**2 + VAR
   80    CONTINUE
         VAR = VAR/DBLE(N)
         XKAPPA = 0.0D0
         XMU2 = XMU*XMU
         IF (XMU2.GT.0.0D0) XKAPPA = 1.0D0 + DBLE(M)/DBLE(N)*VAR/XMU2
         IF (XMU2.GT.0.0D0) SUM2 = SUM2/XMU2/DBLE(N-M)
         IF (XKAPPA.LE.0.0D0 .OR. SUM2.LE.0.0D0) THEN
            IFAIL = 2
            FACT = 1.0D0
         ELSE
            FACT = (XKAPPA*XKAPPA)*SUM2
            FACT = FACT*SIGMA*SIGMA
         END IF
         DO 120 J = 1, M
            J1 = J + 1
            IJ = (J-1)*(M+1) + J1
            C(J,J) = WK(IJ)*FACT
            DO 100 I = J1, M
               IJ = (J-1)*(M+1) + I + 1
               C(I,J) = WK(IJ)*FACT
               C(J,I) = C(I,J)
  100       CONTINUE
  120    CONTINUE
      ELSE
C
C        FOR MALLOWS AND SCHWEPPE CASES
C
         IF (INDC.NE.1) THEN
C
C           IF INDC ne 1
C
            IF (INDW.LT.0) THEN
C
C              MALLOWS CASE
C
               DO 140 I = 1, N
                  Z = RS(I)/SIGMA
                  D(I) = G02HAU(Z,IPSI,CPSI,H1,H2,H3)*WGT(I)
                  E(I) = (G07DBX(Z,IPSI,CPSI,H1,H2,H3)*WGT(I))**2
  140          CONTINUE
            ELSE
C
C              SCHWEPPE CASE
C
               DO 160 I = 1, N
                  IF (WGT(I).LE.0.0D0) THEN
                     D(I) = -1.0D0
                     E(I) = 0.0D0
                  ELSE
                     Z = RS(I)/SIGMA/WGT(I)
                     D(I) = G02HAU(Z,IPSI,CPSI,H1,H2,H3)
                     E(I) = (G07DBX(Z,IPSI,CPSI,H1,H2,H3)*WGT(I))**2
                  END IF
  160          CONTINUE
            END IF
         ELSE
C
C           IF INDC eq 1
C
            IF (INDW.LT.0) THEN
C
C              MALLOWS CASE
C
               S1 = 0.0D0
               S2 = 0.0D0
               DO 180 J = 1, N
                  IF (WGT(J).GT.0.0D0) THEN
                     Z = RS(J)/SIGMA
                     S1 = G02HAU(Z,IPSI,CPSI,H1,H2,H3) + S1
                     S2 = G02HAU(Z,IPSI,CPSI,H1,H2,H3)**2 + S2
                  END IF
  180          CONTINUE
               DO 200 I = 1, N
                  D(I) = S1/DBLE(N)*WGT(I)
                  E(I) = S2/DBLE(N)*WGT(I)*WGT(I)
  200          CONTINUE
            ELSE
C
C              SCHWEPPE CASE
C
               DO 240 I = 1, N
                  S1 = 0.0D0
                  S2 = 0.0D0
                  IF (WGT(I).GT.0.0D0) THEN
                     DO 220 J = 1, N
                        Z = RS(J)/SIGMA/WGT(I)
                        S1 = G02HAU(Z,IPSI,CPSI,H1,H2,H3) + S1
                        S2 = G07DBX(Z,IPSI,CPSI,H1,H2,H3)**2 + S2
  220                CONTINUE
                  END IF
                  D(I) = S1/DBLE(N)
                  E(I) = S2/DBLE(N)*WGT(I)*WGT(I)
  240          CONTINUE
            END IF
         END IF
C
C        CALCULATE ESTIMATE OF THE FORM INV(S1)*S2*INV(S1)
C
         RN = DBLE(N)
         FACT = (SIGMA*SIGMA)/RN
C
C        CALCULATE S1
C
         DO 260 I = 1, M
            IJ = (I-1)*N + 1
            CALL DCOPY(N,X(1,I),1,WK(IJ),1)
            CALL F06FCF(N,D,1,WK(IJ),1)
  260    CONTINUE
         DO 300 J = 1, M
            J1 = J - 1
            IJ = (J-1)*N + 1
            C(J,J) = DDOT(N,X(1,J),1,WK(IJ),1)/RN
            DO 280 I = 1, J1
               C(I,J) = DDOT(N,X(1,I),1,WK(IJ),1)/RN
               C(J,I) = C(I,J)
  280       CONTINUE
  300    CONTINUE
C
C        CALCULATE INV(S1)
C
         CALL F01AAZ(C,IC,M,WK(1),M,WK(MM1),IFAIL2)
         IF (IFAIL2.NE.0) THEN
            IFAIL = 1
            RETURN
         END IF
C
C        CALCULATE S2
C
         DO 320 I = 1, M
            IJ = (I-1)*N + MM1
            CALL DCOPY(N,X(1,I),1,WK(IJ),1)
            CALL F06FCF(N,E,1,WK(IJ),1)
  320    CONTINUE
         DO 360 J = 1, M
            DO 340 I = 1, J
               IJ = (J-1)*N + MM1
               C(I,J) = DDOT(N,X(1,I),1,WK(IJ),1)/RN
  340       CONTINUE
  360    CONTINUE
C
C        CALCULATE INV(S1)*S2
C
         DO 400 J = 1, M
            DO 380 I = 1, M
               IJ = (J-1)*M + I + MM1
               WK(IJ) = DDOT(J,WK(I),M,C(1,J),1)
               IF (J.NE.M) WK(IJ) = WK(IJ) + DDOT(M-J,WK(I+J*M),M,C(J,
     *                              J+1),IC)
  380       CONTINUE
  400    CONTINUE
C
C        CALCULATE INV(S1)*S2*INV(S1)
C
         DO 440 J = 1, M
            J1 = J - 1
            IJ = (J-1)*M + 1
            C(J,J) = DDOT(M,WK(J+MM1),M,WK(IJ),1)*FACT
            DO 420 I = 1, J1
               C(I,J) = DDOT(M,WK(I+MM1),M,WK(IJ),1)*FACT
               C(J,I) = C(I,J)
  420       CONTINUE
  440    CONTINUE
      END IF
      RETURN
C
      END
