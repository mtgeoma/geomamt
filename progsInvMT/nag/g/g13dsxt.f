      SUBROUTINE G13DSX(K,M,P,Q,PAR,NPAR,DEL,X,TEMP,HOLD,ARRAY,TEMPO,IM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IM, K, M, NPAR, P, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  ARRAY(K,M*K), DEL(K,K), HOLD(IM,M*K*K),
     *                  PAR(NPAR), TEMP(M*K*K,M*K*K), TEMPO(K,K),
     *                  X(NPAR,M*K*K)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           A, B, G, H, I, J, J2, K3, L, L2, MK2, PK2
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     This subroutine builds up the X' matrix (stored as X)
C
C     initialise X array to zero
C
      K3 = K*K
      PK2 = P*K3
      MK2 = M*K3
      DO 40 J = 1, MK2
         DO 20 I = 1, NPAR
            X(I,J) = 0.0D0
   20    CONTINUE
   40 CONTINUE
C
C     Set up the autoregressive part of X first
C
C     Initialise ARRAY and HOLD to zero
C
      IF (P.EQ.0) GO TO 600
      DO 80 J = 1, M*K
         DO 60 I = 1, K
            ARRAY(I,J) = 0.0D0
            HOLD(I,J) = 0.0D0
   60    CONTINUE
   80 CONTINUE
C
C                      -1
C     put HOLD = PHI(B)  * THETA(B)
C
C     first set HOLD = THETA(B)
C
      DO 140 L = 1, Q
         DO 120 J = 1, K
            DO 100 I = 1, K
               HOLD(I,(L-1)*K+J) = -PAR(PK2+(L-1)*K3+(I-1)*K+J)
  100       CONTINUE
  120    CONTINUE
  140 CONTINUE
C
      DO 240 L = 1, M
         DO 220 I = 1, K
            DO 200 J = 1, K
               SUM = 0.0D0
               DO 180 L2 = 1, MIN(P,L)
                  IF (L.GT.L2) THEN
                     DO 160 A = 1, K
                        SUM = SUM + PAR((L2-1)*K3+(I-1)*K+A)*HOLD(A,
     *                        (L-L2-1)*K+J)
  160                CONTINUE
                  ELSE
                     SUM = SUM + PAR((L2-1)*K3+(I-1)*K+J)
                  END IF
  180          CONTINUE
               HOLD(I,(L-1)*K+J) = HOLD(I,(L-1)*K+J) + SUM
  200       CONTINUE
  220    CONTINUE
  240 CONTINUE
C
      DO 580 A = 1, K
         DO 560 B = 1, K
C
C           Put ARRAY = D(a,b) * HOLD * DEL
C
            DO 300 L = 1, M
               DO 280 J = 1, K
                  SUM = 0.0D0
                  DO 260 J2 = 1, K
                     SUM = SUM + HOLD(B,(L-1)*K+J2)*DEL(J2,J)
  260             CONTINUE
                  ARRAY(A,(L-1)*K+J) = SUM
  280          CONTINUE
  300       CONTINUE
C
C           set TEMPO = D(a,b) * DEL
C
            DO 340 J = 1, K
               DO 320 I = 1, K
                  TEMPO(I,J) = 0.0D0
  320          CONTINUE
               TEMPO(A,J) = DEL(B,J)
  340       CONTINUE
C
C                                  -1
C           Now put TEMP = THETA(B)   * (TEMPO : ARRAY)
C
            DO 460 L = 1, M
               DO 440 I = 1, K
                  DO 420 J = 1, K
                     SUM = ARRAY(I,(L-1)*K+J)
                     DO 400 L2 = 1, MIN(Q,L)
                        IF (L.GT.L2) THEN
                           DO 360 H = 1, K
                              SUM = SUM + PAR(PK2+(L2-1)*K3+(I-1)*K+H)
     *                              *TEMP(H,(L-L2-1)*K+J)
  360                      CONTINUE
                        ELSE
                           DO 380 H = 1, K
                              SUM = SUM + PAR(PK2+(L2-1)*K3+(I-1)*K+H)
     *                              *TEMPO(H,J)
  380                      CONTINUE
                        END IF
  400                CONTINUE
                     TEMP(I,(L-1)*K+J) = SUM
  420             CONTINUE
  440          CONTINUE
  460       CONTINUE
C
C           Assign to X array
C
            DO 540 G = 1, K
               DO 520 H = 1, K
                  DO 500 I = 1, P
                     X((I-1)*K3+(A-1)*K+B,(I-1)*K3+(G-1)*K+H) = TEMPO(G,
     *                 H)
                     DO 480 J = I + 1, M
                        X((I-1)*K3+(A-1)*K+B,(J-1)*K3+(G-1)*K+H)
     *                    = TEMP(G,(J-I-1)*K+H)
  480                CONTINUE
  500             CONTINUE
  520          CONTINUE
  540       CONTINUE
  560    CONTINUE
  580 CONTINUE
C
C
C     Now set up the moving average part of X
C
C     Initialise HOLD to zero
C
  600 IF (Q.EQ.0) RETURN
      DO 640 J = 1, M*K
         DO 620 I = 1, K
            HOLD(I,J) = 0.0D0
  620    CONTINUE
  640 CONTINUE
C                        -1
C     put HOLD = THETA(B)
C
      DO 740 L = 1, M
         DO 720 I = 1, K
            DO 700 J = 1, K
               SUM = 0.0D0
               DO 680 L2 = 1, MIN(Q,L)
                  IF (L.GT.L2) THEN
                     DO 660 A = 1, K
                        SUM = SUM + PAR(PK2+(L2-1)*K3+(I-1)*K+A)*HOLD(A,
     *                        (L-L2-1)*K+J)
  660                CONTINUE
                  ELSE
                     SUM = SUM + PAR(PK2+(L2-1)*K3+(I-1)*K+J)
                  END IF
  680          CONTINUE
               HOLD(I,(L-1)*K+J) = SUM
  700       CONTINUE
  720    CONTINUE
  740 CONTINUE
C
      DO 960 A = 1, K
         DO 940 B = 1, K
C
C           put TEMP = HOLD * D(a,b) * DEL
C
            DO 800 L = 1, M
               DO 780 J = 1, K
                  SUM = DEL(B,J)
                  DO 760 I = 1, K
                     TEMP(I,(L-1)*K+J) = HOLD(I,(L-1)*K+A)*SUM
  760             CONTINUE
  780          CONTINUE
  800       CONTINUE
C
C           set TEMPO = D(a,b) * DEL
C
            DO 840 J = 1, K
               DO 820 I = 1, K
                  TEMPO(I,J) = 0.0D0
  820          CONTINUE
               TEMPO(A,J) = DEL(B,J)
  840       CONTINUE
C
C           Assign to X array
C
            DO 920 G = 1, K
               DO 900 H = 1, K
                  DO 880 I = 1, Q
                     X(PK2+(I-1)*K3+(A-1)*K+B,(I-1)*K3+(G-1)*K+H)
     *                 = TEMPO(G,H)
                     DO 860 J = I + 1, M
                        X(PK2+(I-1)*K3+(A-1)*K+B,(J-1)*K3+(G-1)*K+H)
     *                    = TEMP(G,(J-I-1)*K+H)
  860                CONTINUE
  880             CONTINUE
  900          CONTINUE
  920       CONTINUE
  940    CONTINUE
  960 CONTINUE
      RETURN
C
      END
