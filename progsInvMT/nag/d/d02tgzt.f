      SUBROUTINE D02TGZ(N,M,L,X0,X1,K1,KP,B,IB,COEFFX,BDYCX,COEFF,BDYC,
     *                  CF,BC,A,N1,IA1,RB,IRB,IW,J)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-312 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     BC, BDYC, BDYCX, COEFF, COEFFX
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X0, X1
      INTEGER           IA1, IB, IRB, J, K1, KP, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N1,IA1), B(IB,N), RB(IRB,7)
      INTEGER           IW(IRB), L(N), M(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF
      EXTERNAL          CF
C     .. Subroutine Arguments ..
      EXTERNAL          BC, BDYC, BDYCX, COEFF, COEFFX
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, G, PW, RHS, T, WS, X, Y, Z
      INTEGER           BW, I, IM, IQ, IRBJ1, IRBV, J1, J2, M1, MJ, MP,
     *                  MU, P, Q, U, V, W
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF
      EXTERNAL          X01AAF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          F04AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, COS, INT
C     .. Executable Statements ..
      J = 1
      M1 = 0
      DO 20 I = 1, N
         IF (M(I).LE.0) GO TO 580
         M1 = MAX(M1,M(I))
   20 CONTINUE
      M1 = M1 + 1
      Q = N*K1
      IF (IB.LT.K1 .OR. K1.LT.M1 .OR. N1.LT.Q .OR. X1.LE.X0) GO TO 580
      WS = KP - 1
      EPS = X02AJF()
      PW = X01AAF(0.0D0)/WS
      G = 2.0D0/(X1-X0)
      BW = N*KP
      DO 280 I = 1, N
         MP = KP + L(I)
         DO 260 P = 1, MP
            RHS = 0.0D0
            DO 60 V = 1, M1
               IRBV = IRB + V + 1
               DO 40 U = 1, N
                  A(U,IRBV) = 0.0D0
   40          CONTINUE
   60       CONTINUE
            IF (P.GT.KP) GO TO 80
            WS = P - 1
            T = MAX(-1.0D0,MIN(1.0D0,COS(PW*WS)))
            X = T/G + (X0+X1)/2.0D0
            CALL COEFFX(X,I,A(1,IRB+2),N1,IA1-IRB-1,RHS,COEFF,CF)
            W = (I-1)*KP + P
            GO TO 100
   80       V = P - KP
            CALL BDYCX(X,I,V,A(1,IRB+2),N1,IA1-IRB-1,RHS,BDYC,BC)
            T = G*0.5D0*((X-X0)+(X-X1))
            BW = BW + 1
            W = BW
  100       RB(1,2) = 1.0D0
            RB(2,2) = T
            IF (K1.LT.3) GO TO 140
            DO 120 J = 3, K1
               RB(J,2) = 2.0D0*T*RB(J-1,2) - RB(J-2,2)
  120       CONTINUE
  140       RB(1,2) = 0.5D0
            DO 240 U = 1, N
               MU = M(U)
               X = 1.0D0
               T = A(U,IRB+2)
               V = (U-1)*K1
               DO 160 Q = 1, K1
                  RB(Q,3) = RB(Q,2)*T
  160          CONTINUE
               DO 200 J = 1, MU
                  J1 = J + V
                  A(W,J1) = RB(1,3)
                  X = X*G
                  IRBJ1 = IRB + J + 2
                  T = A(U,IRBJ1)*X
                  Y = 0.0D0
                  Z = 0.0D0
                  J1 = K1 - J
                  DO 180 Q = 1, J1
                     WS = 2*Q
                     RB(Q,3) = RB(Q+1,3)/WS - Y + RB(Q,2)*T
                     Y = Z
                     Z = RB(Q+1,3)/WS
  180             CONTINUE
  200          CONTINUE
               MJ = MU + 1
               DO 220 J = MJ, K1
                  J1 = J + V
                  J2 = J - MU
                  A(W,J1) = RB(J2,3)
  220          CONTINUE
  240       CONTINUE
            A(W,IRB+1) = RHS
  260    CONTINUE
  280 CONTINUE
      V = N*KP
      Q = N*K1
      I = N1
  300 G = 0.0D0
      DO 320 J = 1, Q
         IF (ABS(A(I,J)).LE.G) GO TO 320
         G = ABS(A(I,J))
         P = J
  320 CONTINUE
      J = 3
      IF (G.LT.EPS**0.75D0) GO TO 580
      DO 340 J = 1, I
         T = A(J,P)
         A(J,P) = A(J,Q)
         A(J,Q) = T
  340 CONTINUE
      IQ = Q - 1
      DO 360 J = 1, IQ
         A(I,J) = A(I,J)/A(I,Q)
  360 CONTINUE
      A(I,IRB+1) = A(I,IRB+1)/A(I,Q)
      A(I,Q) = P
      IM = I - 1
      DO 400 P = 1, IM
         T = A(P,Q)
         DO 380 J = 1, IQ
            A(P,J) = A(P,J) - T*A(I,J)
  380    CONTINUE
         A(P,IRB+1) = A(P,IRB+1) - T*A(I,IRB+1)
  400 CONTINUE
      Q = Q - 1
      I = I - 1
      IF (I.GT.V) GO TO 300
      J = 1
      CALL F04AMF(A,N1,RB,IRB,A(1,IRB+1),N1,V,Q,1,EPS,A(1,IRB+2)
     *            ,N1,RB(1,4),RB(1,5),RB(1,6),RB(1,7),A(1,2*IRB+2),IW,J)
      IF (J.EQ.0) GO TO 420
      J = J + 2
      GO TO 580
  420 IM = V + 1
      DO 460 I = IM, N1
         T = A(I,IRB+1)
         DO 440 J = 1, Q
            T = T - A(I,J)*RB(J,1)
  440    CONTINUE
         Q = Q + 1
         P = INT(A(I,Q)+0.01D0)
         RB(Q,1) = RB(P,1)
         RB(P,1) = T
  460 CONTINUE
      DO 560 I = 1, N
         U = (I-1)*K1
         MU = M(I)
         J1 = K1 - MU
         DO 480 J = 1, J1
            IM = U + MU + J
            A(J,MU+1) = RB(IM,1)
  480    CONTINUE
         J = MU
  500    V = K1 - J + 1
  520    T = 0.0D0
         IF (V.LT.K1-J) T = A(V+1,J+1)
         WS = 2*V - 2
         A(V,J) = (A(V-1,J+1)-T)/WS
         V = V - 1
         IF (V.GT.1) GO TO 520
         IM = U + J
         A(1,J) = RB(IM,1)
         J = J - 1
         IF (J.GT.0) GO TO 500
         DO 540 J = 1, K1
            B(J,I) = A(J,1)
  540    CONTINUE
  560 CONTINUE
      J = 0
  580 RETURN
      END
