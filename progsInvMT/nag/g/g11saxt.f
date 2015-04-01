      SUBROUTINE G11SAX(MAXIT,IPRINT,A,C,XN,AX,Q,N2,S,X,RL,LP,GPROB,
     *                  GRATOL,NOMON,ROOTPI,N3,IERROR,R,G,P,PHI,NROWXR,
     *                  NITER)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CARRY OUT UP TO MAXIT ITERATIONS OF THE E-M ALGORITHM
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GRATOL, ROOTPI
      INTEGER           IERROR, IPRINT, LP, MAXIT, N2, N3, NITER,
     *                  NROWXR, Q, S
      LOGICAL           GPROB, NOMON
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N2), AX(20), C(N2), G(N3), P(S), PHI(N2,20),
     *                  R(N2,20), XN(20)
      INTEGER           RL(S)
      LOGICAL           X(NROWXR,N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A2, B2, GPJNRM, LIMIT, P2, P7, PROD, S1, S2, S3,
     *                  S4, SSM, SUM, SUM2, TOLL, TT, X2, XBAR, YBAR
      INTEGER           I, I2, IFAIL, ITN, J, K, L
      LOGICAL           CGE
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  LX(20), NK(20), W(20), Y(20), ZX(20)
C     .. External Functions ..
      DOUBLE PRECISION  G01CEF, S15ABF, X02AJF, X02AMF
      EXTERNAL          G01CEF, S15ABF, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          G11SAS, G11SAW, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, DBLE
C     .. Executable Statements ..
      LIMIT = LOG(X02AMF())
      IERROR = 0
      S3 = X02AJF()
      IFAIL = 0
      TOLL = G01CEF(S3,IFAIL)
C
      DO 520 ITN = 1, MAXIT
C
         NITER = ITN
C
C        THE EXPECTATION STEP
C
C        CALCULATE THE ELEMENTS OF PHI OVER ALL ITEMS AND NODES
C
         IF (GPROB) THEN
            DO 40 J = 1, N2
               DO 20 K = 1, Q
                  IFAIL = 0
                  SSM = A(J)*XN(K) + C(J)
                  PHI(J,K) = S15ABF(SSM,IFAIL)
                  R(J,K) = 0.0D0
   20          CONTINUE
   40       CONTINUE
C
         ELSE
            DO 120 J = 1, N2
               DO 100 K = 1, Q
                  SUM = A(J)*XN(K) + C(J)
                  R(J,K) = 0.0D0
                  IF (SUM.LE.0.0D0) THEN
                     IF (SUM.LT.LIMIT) THEN
                        SUM2 = 0
                        GO TO 60
                     END IF
                     SUM2 = EXP(SUM)
   60                PHI(J,K) = SUM2/(1.0D0+SUM2)
                  ELSE
                     IF ((-SUM).LT.(LIMIT)) THEN
                        PHI(J,K) = 1.0D0
                        GO TO 80
                     END IF
                     PHI(J,K) = 1.0D0/(1.0D0+EXP(-SUM))
   80             END IF
C
  100          CONTINUE
  120       CONTINUE
C
         END IF
C
         DO 140 K = 1, Q
            NK(K) = 0.0D0
  140    CONTINUE
C
C        CALCULATE R(J,K) AND NK(K) FOR J=1,2,,,N2   ,   K=1,2,,,Q
C
         DO 260 L = 1, S
C
            P7 = 0.0D0
            DO 180 K = 1, Q
C
C              CALCULATE LX(L,K)
C
               PROD = 1.0D0
C
               DO 160 J = 1, N2
C
                  IF (X(L,J)) THEN
                     PROD = PROD*PHI(J,K)
                  ELSE
                     PROD = PROD*(1.0D0-PHI(J,K))
                  END IF
C
  160          CONTINUE
               LX(K) = PROD
               P7 = P7 + PROD*AX(K)
  180       CONTINUE
C
C           CALCULATE P(L)
C
            P(L) = 1.0D0/P7
            P7 = 1.0D0/P7
C
C           UPDATE R AND NK
C
            DO 220 J = 1, N2
               DO 200 K = 1, Q
                  IF (X(L,J)) R(J,K) = R(J,K) + DBLE(RL(L))*LX(K)*P7
  200          CONTINUE
  220       CONTINUE
C
            DO 240 K = 1, Q
               NK(K) = NK(K) + DBLE(RL(L))*LX(K)*P7
  240       CONTINUE
C
  260    CONTINUE
C
C        CORRECT ELEMENTS OF R AND NK
C
         DO 300 J = 1, N2
            DO 280 K = 1, Q
               R(J,K) = R(J,K)*AX(K)
  280       CONTINUE
  300    CONTINUE
C
         DO 320 K = 1, Q
            NK(K) = NK(K)*AX(K)
  320    CONTINUE
C
C        CARRY OUT THE MAXIMISATION STEP
C
         DO 440 J = 1, N2
C
C           CALCULATE WEIGHTS AND (WORKING PROBITS * WEIGHTS)
C
            DO 400 K = 1, Q
C
               ZX(K) = A(J)*XN(K) + C(J)
               SUM = ZX(K)
               IF (GPROB) THEN
                  IF ((-(SUM*SUM)/2.0D0).LT.(LIMIT)) THEN
                     SUM2 = 0.0D0
                     GO TO 340
                  END IF
                  SUM2 = ROOTPI*EXP(-(SUM*SUM)/2.0D0)
  340             IF (SUM.GE.-TOLL) THEN
C
                     S3 = SUM*SUM
                     S3 = SUM*((S3-0.5D0)/(S3-1.0D0))
                     P7 = (SUM2/PHI(J,K))*S3
                     W(K) = P7
                     S4 = (R(J,K)*S3)/(PHI(J,K)*NK(K))
                     Y(K) = P7*SUM + S4 - S3
C
                  ELSE
                     IF (SUM.LE.TOLL) THEN
C
                        S3 = SUM*SUM
                        S3 = SUM*((S3-0.5D0)/(S3-1.0D0))
                        P7 = -S3*SUM2/(1.0D0-PHI(J,K))
                        W(K) = P7
                        S4 = (R(J,K)/NK(K)-PHI(J,K))*S3
                        Y(K) = P7*SUM - (S4/(1.0D0-PHI(J,K)))
C
                     ELSE
C
                        SSM = SUM
                        P7 = SUM2/PHI(J,K)
                        IFAIL = 0
                        TT = S15ABF(-SSM,IFAIL)
                        PROD = SUM2/TT
                        W(K) = P7*PROD
                        P2 = W(K)*SUM - (SUM2/TT)
                        Y(K) = P2 + ((R(J,K)/NK(K))*((SUM2/PHI(J,K))/TT)
     *                         )
C
                     END IF
C
                  END IF
C
               ELSE
C
                  IF (SUM.LE.0.0D0) THEN
                     IF ((SUM).LT.(LIMIT)) THEN
                        SUM2 = 0
                        GO TO 360
                     END IF
                     SUM2 = EXP(SUM)
  360                W(K) = SUM2/((1.0D0+SUM2)**2)
                  ELSE
                     IF ((-SUM).LT.(LIMIT)) THEN
                        SUM2 = 0
                        GO TO 380
                     END IF
                     SUM2 = EXP(-SUM)
  380                W(K) = SUM2/((1.0D0+SUM2)**2)
                  END IF
                  Y(K) = W(K)*SUM + (R(J,K)/NK(K)-PHI(J,K))
C
               END IF
  400       CONTINUE
C
C           CALCULATE SNW,SNWX,SNWY,,, ETC.
C
            B2 = 0.0D0
            SUM = 0.0D0
            X2 = 0.0D0
            A2 = 0.0D0
            PROD = 0.0D0
            S1 = 0.0D0
            S2 = 0.0D0
C
            DO 420 K = 1, Q
C
               SUM2 = NK(K)*W(K)
               XBAR = XN(K)*SUM2
               YBAR = XBAR*XN(K)
               B2 = B2 + XBAR
               A2 = A2 + SUM2
               PROD = PROD + NK(K)*Y(K)
               SUM = SUM + NK(K)*XN(K)*Y(K)
               X2 = X2 + YBAR
               SUM2 = NK(K)*(Y(K)-ZX(K)*W(K))
               S1 = S1 + SUM2
               S2 = S2 + SUM2*XN(K)
C
  420       CONTINUE
C
            XBAR = B2/A2
            YBAR = PROD/A2
C
C           UPDATE A(J) AND C(J) AND THEN TEST FOR CONVERGENCE
C
            SUM = SUM - XBAR*PROD
            A(J) = SUM/(X2-XBAR*XBAR*A2)
            C(J) = YBAR - A(J)*XBAR
C
            G(2*J-1) = S2
            G(2*J) = S1
C
  440    CONTINUE
C
         GPJNRM = ABS(G(1))
         DO 460 I = 2, N3
            IF (ABS(G(I)).GT.GPJNRM) GPJNRM = ABS(G(I))
  460    CONTINUE
C
         CGE = GPJNRM .LT. GRATOL
         IF (NOMON .OR. (IPRINT.EQ.0)) GO TO 480
         I2 = ITN/IPRINT*IPRINT
         IF ((I2.EQ.ITN) .OR. CGE) CALL G11SAS(N2,A,C,LP,P,RL,S,G,ITN,
     *                                  GPJNRM,N3)
C
  480    IF (CGE .AND. (Q.EQ.20)) THEN
            IF (IPRINT.EQ.0) CALL G11SAS(N2,A,C,LP,P,RL,S,G,ITN,GPJNRM,
     *                                   N3)
            RETURN
         END IF
         IF (CGE .AND. (Q.EQ.10)) THEN
            Q = 20
         END IF
         IF (CGE) THEN
            IF (IPRINT.GT.0) THEN
               WRITE (REC,FMT=99999) Q
               CALL X04BAF(LP,REC)
            END IF
            CALL G11SAW(Q,XN,AX,LP,IPRINT,ROOTPI)
         END IF
C
         DO 500 I = 1, N2
            IF (ABS(A(I)).GE.10.0D0) THEN
               IERROR = 5
               Q = 20
               CALL G11SAW(Q,XN,AX,LP,IPRINT,ROOTPI)
               RETURN
            END IF
  500    CONTINUE
C
  520 CONTINUE
C
      IERROR = 3
C
      RETURN
C
99999 FORMAT (' THE NUMBER OF QUADRATURE POINTS HAS NOW BEEN INCREASED',
     *  ' TO ',I2)
      END
