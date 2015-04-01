      SUBROUTINE G11SAZ(A,C,GPROB,N2,LL,VAR,G,XN,AX,Q,S,X,RL,ALPHA,
     *                  GAMMA,RR,ROOTPI,N3,PHI,P,NROWXR,IAA,PL,R2,M,IAN,
     *                  IERROR)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATE VARIANCE-COVARIANCE MATRIX OF ITEM PARAMETER
C     ESTIMATES
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  LL, ROOTPI
      INTEGER           IAA, IAN, IERROR, N2, N3, NROWXR, Q, S
      LOGICAL           GPROB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N2), ALPHA(N2), AX(20), C(N2), G(N3),
     *                  GAMMA(N2), M(IAN,N3), P(S), PHI(N2,20), PL(N3),
     *                  R2(N2,20), RR(N2,20), VAR(IAA,N3), XN(20)
      INTEGER           RL(S)
      LOGICAL           X(NROWXR,N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A2, A3, A7, LIMIT, PROD, SUM, SUM2, SUM3, SUM4,
     *                  SUM5
      INTEGER           I, IFAIL, J, J4, K, L
C     .. Local Arrays ..
      DOUBLE PRECISION  LX(20)
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X02AMF
      EXTERNAL          S15ABF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          F01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, SQRT, DBLE
C     .. Executable Statements ..
      LIMIT = LOG(X02AMF())
C
C     CALCULATE THE ELEMENTS OF PHI USING A STANDARD ROUTINE
C
      IF (GPROB) THEN
         DO 60 J = 1, N2
            DO 40 K = 1, Q
               A2 = C(J) + A(J)*XN(K)
               IFAIL = 0
               PHI(J,K) = S15ABF(A2,IFAIL)
               IF ((-(A2*A2)/2.0D0).LT.(LIMIT)) THEN
                  RR(J,K) = 0
                  GO TO 20
               END IF
               RR(J,K) = ROOTPI*EXP(-(A2*A2)/2.0D0)
   20          R2(J,K) = -A2
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 140 J = 1, N2
            DO 120 K = 1, Q
               SUM = C(J) + A(J)*XN(K)
               IF (SUM.LE.0.0D0) THEN
                  IF ((SUM).LT.(LIMIT)) THEN
                     SUM2 = 0
                     GO TO 80
                  END IF
                  SUM2 = EXP(SUM)
   80             PHI(J,K) = SUM2/(1.0D0+SUM2)
                  RR(J,K) = SUM2/((1.0D0+SUM2)**2)
                  R2(J,K) = (1.0D0-SUM2)/(1.0D0+SUM2)
               ELSE
                  IF ((-SUM).LT.(LIMIT)) THEN
                     SUM2 = 0
                     GO TO 100
                  END IF
                  SUM2 = EXP(-SUM)
  100             PHI(J,K) = 1.0D0/(1.0D0+SUM2)
                  RR(J,K) = SUM2/((1.0D0+SUM2)**2)
                  R2(J,K) = (SUM2-1.0D0)/(1.0D0+SUM2)
               END IF
  120       CONTINUE
  140    CONTINUE
C
      END IF
C
C     INITIALISE M AND G TO ZERO
C
      DO 180 I = 1, N3
         DO 160 J = 1, N3
            M(I,J) = 0.0D0
  160    CONTINUE
         G(I) = 0.0D0
  180 CONTINUE
C
      LL = 0.0D0
C
      DO 600 L = 1, S
C
C        CALCULATE LX(K) AND P(L)
C
         DO 220 K = 1, Q
            PROD = 1.0D0
            DO 200 J = 1, N2
C
               IF (X(L,J)) THEN
                  PROD = PROD*PHI(J,K)
               ELSE
                  PROD = PROD*(1.0D0-PHI(J,K))
               END IF
C
  200       CONTINUE
            LX(K) = PROD
  220    CONTINUE
C
C        CALCULATE P(L)
C
         SUM = 0.0D0
         DO 240 K = 1, Q
            SUM = SUM + LX(K)*AX(K)
  240    CONTINUE
         P(L) = SUM
C
C
         DO 340 J = 1, N2
            SUM = 0.0D0
            SUM2 = 0.0D0
            SUM3 = 0.0D0
            SUM4 = 0.0D0
            SUM5 = 0.0D0
            IF (X(L,J)) THEN
               DO 280 K = 1, Q
                  IF (PHI(J,K).LE.0.0D0) THEN
                     A3 = AX(K)*RR(J,K)
                     DO 260 J4 = 1, N2
                        IF (J4.NE.J) THEN
                           IF (X(L,J4)) THEN
                              A3 = A3*PHI(J4,K)
                           ELSE
                              A3 = A3*(1.0D0-PHI(J4,K))
                           END IF
                        END IF
  260                CONTINUE
C
                  ELSE
                     A3 = (AX(K)/PHI(J,K))*LX(K)*RR(J,K)
                  END IF
                  SUM = SUM + A3
                  A7 = A3*XN(K)
                  SUM2 = SUM2 + A7
                  SUM3 = SUM3 + A7*R2(J,K)*XN(K)
                  SUM4 = SUM4 + A3*R2(J,K)
                  SUM5 = SUM5 + A7*R2(J,K)
  280          CONTINUE
C
            ELSE
C
               DO 320 K = 1, Q
                  IF (PHI(J,K).GE.1.0D0) THEN
                     A3 = -AX(K)*RR(J,K)
                     DO 300 J4 = 1, N2
                        IF (J4.NE.J) THEN
                           IF (X(L,J4)) THEN
                              A3 = A3*PHI(J4,K)
                           ELSE
                              A3 = A3*(1.0D0-PHI(J4,K))
                           END IF
                        END IF
  300                CONTINUE
                  ELSE
                     A3 = -(AX(K)/(1.0D0-PHI(J,K)))*LX(K)*RR(J,K)
                  END IF
C
                  SUM = SUM + A3
                  A7 = A3*XN(K)
                  SUM2 = SUM2 + A7
                  SUM3 = SUM3 + A7*R2(J,K)*XN(K)
                  SUM4 = SUM4 + A3*R2(J,K)
                  SUM5 = SUM5 + A7*R2(J,K)
C
  320          CONTINUE
C
            END IF
C
            PL(2*J-1) = SUM2
            PL(2*J) = SUM
            VAR(2*J-1,2*J-1) = SUM3
            VAR(2*J,2*J) = SUM4
            VAR(2*J-1,2*J) = SUM5
  340    CONTINUE
C
C        CALCULATE NON-DIAGONAL ELEMENTS OF 2ND DERIVATIVES OF P(L)'S
C        AND STORE AS VAR
C
         DO 540 I = 1, N2 - 1
            DO 520 J = I + 1, N2
               SUM = 0.0D0
               SUM2 = 0.0D0
               SUM3 = 0.0D0
C
               IF (X(L,J)) THEN
                  IF (X(L,I)) THEN
                     DO 380 K = 1, Q
                        IF ((PHI(J,K).LE.0.0D0) .OR. (PHI(I,K).LE.0.0D0)
     *                      ) THEN
                           A3 = AX(K)*RR(J,K)*RR(I,K)
                           DO 360 J4 = 1, N2
                              IF ((J4.NE.I) .AND. (J4.NE.J)) THEN
                                 IF (X(L,J4)) THEN
                                    A3 = A3*PHI(J4,K)
                                 ELSE
                                    A3 = A3*(1.0D0-PHI(J4,K))
                                 END IF
                              END IF
  360                      CONTINUE
                        ELSE
                           A3 = AX(K)*(RR(J,K)/PHI(J,K))*(RR(I,K)/PHI(I,
     *                          K))*LX(K)
                        END IF
C
                        SUM = SUM + A3
                        SUM2 = SUM2 + A3*XN(K)
                        SUM3 = SUM3 + A3*XN(K)*XN(K)
C
  380                CONTINUE
                  ELSE
                     DO 420 K = 1, Q
                        IF ((PHI(J,K).LE.0.0D0) .OR. (PHI(I,K).GE.1.0D0)
     *                      ) THEN
                           A3 = -AX(K)*RR(J,K)*RR(I,K)
                           DO 400 J4 = 1, N2
                              IF ((J4.NE.I) .AND. (J4.NE.J)) THEN
                                 IF (X(L,J4)) THEN
                                    A3 = A3*PHI(J4,K)
                                 ELSE
                                    A3 = A3*(1.0D0-PHI(J4,K))
                                 END IF
                              END IF
  400                      CONTINUE
C
                        ELSE
C
                           A3 = -AX(K)*(RR(J,K)/PHI(J,K))*(RR(I,K)
     *                          /(1.0D0-PHI(I,K)))*LX(K)
                        END IF
C
                        SUM = SUM + A3
                        SUM2 = SUM2 + A3*XN(K)
                        SUM3 = SUM3 + A3*XN(K)*XN(K)
  420                CONTINUE
C
                  END IF
C
               ELSE
C
                  IF (X(L,I)) THEN
                     DO 460 K = 1, Q
                        IF ((PHI(J,K).GE.1.0D0) .OR. (PHI(I,K).LE.0.0D0)
     *                      ) THEN
                           A3 = -AX(K)*RR(J,K)*RR(I,K)
                           DO 440 J4 = 1, N2
                              IF ((J4.NE.I) .AND. (J4.NE.J)) THEN
                                 IF (X(L,J4)) THEN
                                    A3 = A3*PHI(J4,K)
                                 ELSE
                                    A3 = A3*(1.0D0-PHI(J4,K))
                                 END IF
                              END IF
  440                      CONTINUE
                        ELSE
                           A3 = -AX(K)*(RR(J,K)/(1.0D0-PHI(J,K)))*LX(K)
     *                          *(RR(I,K)/PHI(I,K))
                        END IF
C
                        SUM = SUM + A3
                        SUM2 = SUM2 + A3*XN(K)
                        SUM3 = SUM3 + A3*XN(K)*XN(K)
C
  460                CONTINUE
                  ELSE
                     DO 500 K = 1, Q
                        IF ((PHI(J,K).GE.1.0D0) .OR. (PHI(I,K).GE.1.0D0)
     *                      ) THEN
                           A3 = AX(K)*RR(J,K)*RR(I,K)
                           DO 480 J4 = 1, N2
                              IF ((J4.NE.I) .AND. (J4.NE.J)) THEN
                                 IF (X(L,J4)) THEN
                                    A3 = A3*PHI(J4,K)
                                 ELSE
                                    A3 = A3*(1.0D0-PHI(J4,K))
                                 END IF
                              END IF
  480                      CONTINUE
C
                        ELSE
C
                           A3 = AX(K)*(RR(J,K)/(1.0D0-PHI(J,K)))*LX(K)
     *                          *(RR(I,K)/(1.0D0-PHI(I,K)))
                        END IF
C
                        SUM = SUM + A3
                        SUM2 = SUM2 + A3*XN(K)
                        SUM3 = SUM3 + A3*XN(K)*XN(K)
  500                CONTINUE
C
                  END IF
C
               END IF
C
               VAR(2*I-1,2*J-1) = SUM3
               VAR(2*I-1,2*J) = SUM2
               VAR(2*I,2*J-1) = SUM2
               VAR(2*I,2*J) = SUM
  520       CONTINUE
  540    CONTINUE
C
C        UPDATE - 2 ND DERIVATIVE MATRIX (UPPER TRIANGLE OF M)
C        AND THE GRADIENT VECTOR
C
         DO 580 J = 1, N3
            DO 560 K = J, N3
               SUM = (PL(J)/P(L))*PL(K)
               SUM = (SUM-VAR(J,K))/P(L)
               M(J,K) = M(J,K) + (DBLE(RL(L))*SUM)
  560       CONTINUE
            G(J) = G(J) + DBLE(RL(L))*(PL(J)/P(L))
  580    CONTINUE
C
C        CALCULATE THE LOG LIKELIHOOD KERNEL
C
         LL = LL + DBLE(RL(L))*LOG(P(L))
C
  600 CONTINUE
C
C     CORRECT COMPONENTS OF GRADIENT VECTOR
C
      IF (GPROB) THEN
         DO 620 J = 1, N2
            SUM = SQRT(1.0D0-ALPHA(J)*ALPHA(J))
            G(2*J-1) = G(2*J-1)/(SUM**3)
            G(2*J) = -G(2*J)/SUM
  620    CONTINUE
      ELSE
         DO 640 J = 1, N2
            G(2*J) = G(2*J)/(GAMMA(J)*(1.0D0-GAMMA(J)))
  640    CONTINUE
      END IF
C
C     CALL INVERSION ROUTINE TO CALCULATE COVARIANCE MATRIX OF
C     A(J)'S AND C(J)'S
C
      IFAIL = 1
      CALL F01ABF(M,IAN,N3,VAR,IAA,PL,IFAIL)
C
      IF (IFAIL.NE.0) THEN
         IF (IERROR.EQ.3) THEN
            IERROR = 6
         ELSE
            IERROR = 7
         END IF
         RETURN
      END IF
C
      DO 680 I = 1, N3
         DO 660 J = 1, I
            M(I,J) = VAR(I,J)
            M(J,I) = M(I,J)
  660    CONTINUE
  680 CONTINUE
C
C     CORRECT COVARIANCE MATRIX
C     AND SET UP LOWER TRIANGLE OF CORRECTED COVARIANCE MATRIX
C
      IF (GPROB) THEN
         DO 700 I = 1, N2
            PL(I) = SQRT(1.0D0+A(I)*A(I))
  700    CONTINUE
C
         DO 740 J = 1, N2
            DO 720 I = J, N2
               IF (I.EQ.J) THEN
                  VAR(2*J-1,2*J-1) = M(2*J-1,2*J-1)/(PL(J)**6)
                  SUM = M(2*J,2*J)/(PL(J)**2)
                  SUM2 = 2.0D0*A(J)*C(J)*M(2*J-1,2*J)/(PL(J)**4)
                  PROD = ((A(J)*C(J))**2)*M(2*J-1,2*J-1)/(PL(J)**6)
                  VAR(2*J,2*J) = SUM - SUM2 + PROD
                  SUM = A(J)*C(J)*M(2*J-1,2*J-1)/(PL(J)**6)
                  SUM2 = M(2*J-1,2*J)/(PL(J)**4)
                  VAR(2*J,2*J-1) = SUM - SUM2
               ELSE
                  VAR(2*I-1,2*J-1) = M(2*I-1,2*J-1)/((PL(I)**3)*(PL(J)
     *                               **3))
                  SUM = M(2*I,2*J)/(PL(I)*PL(J))
                  SUM2 = A(J)*C(J)*M(2*I,2*J-1)/(PL(I)*(PL(J)**3))
                  PROD = A(I)*C(I)*M(2*I-1,2*J)/((PL(I)**3)*PL(J))
                  A2 = A(I)*C(I)*A(J)*C(J)*M(2*I-1,2*J-1)/((PL(I)**3)
     *                 *(PL(J)**3))
                  VAR(2*I,2*J) = SUM - SUM2 - PROD + A2
                  SUM = M(2*I-1,2*J)/((PL(I)**3)*PL(J))
                  PROD = A(J)*C(J)*M(2*I-1,2*J-1)/((PL(I)**3)*(PL(J)**3)
     *                   )
                  VAR(2*I-1,2*J) = PROD - SUM
                  SUM = M(2*J-1,2*I)/((PL(J)**3)*PL(I))
                  PROD = A(I)*C(I)*M(2*J-1,2*I-1)/((PL(J)**3)*(PL(I)**3)
     *                   )
                  VAR(2*I,2*J-1) = PROD - SUM
               END IF
  720       CONTINUE
  740    CONTINUE
C
         RETURN
C
      ELSE
         DO 800 J = 1, N2
            IF (C(J).LE.0.0D0) THEN
               IF ((C(J)).LT.(LIMIT)) THEN
                  SUM = 0
                  GO TO 760
               END IF
               SUM = EXP(C(J))
  760          PL(J) = SUM/((1.0D0+SUM)**2)
            ELSE
               IF ((-C(J)).LT.(LIMIT)) THEN
                  SUM = 0
                  GO TO 780
               END IF
               SUM = EXP(-C(J))
  780          PL(J) = SUM/((1.0D0+SUM)**2)
            END IF
  800    CONTINUE
         DO 840 J = 1, N2
            DO 820 I = J, N2
               IF (I.EQ.J) THEN
                  VAR(2*J-1,2*J-1) = M(2*J-1,2*J-1)
                  VAR(2*J,2*J) = PL(J)*PL(J)*M(2*J,2*J)
                  VAR(2*J,2*J-1) = PL(J)*M(2*J-1,2*J)
               ELSE
                  VAR(2*I-1,2*J-1) = M(2*I-1,2*J-1)
                  VAR(2*I,2*J) = PL(I)*PL(J)*M(2*I,2*J)
                  VAR(2*I-1,2*J) = PL(J)*M(2*I-1,2*J)
                  VAR(2*I,2*J-1) = PL(I)*M(2*J-1,2*I)
               END IF
  820       CONTINUE
  840    CONTINUE
C
         RETURN
      END IF
      END
