      SUBROUTINE D02RAR(M,N,IRN,NIRN,IP,NIP,H,X,F,HMAX,UMIN,UAIM,UMAX,
     *                  EPS,EPS1,ADJUST,A,IG,NIG,W,NW,Z,IND,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     IG(J) POINTS TO THE FIRST ENTRY OF GROUP J
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, EPS1, UAIM, UMAX, UMIN
      INTEGER           IFAIL, IND, M, N, NIG, NIP, NIRN, NW
      LOGICAL           ADJUST
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NIRN), F(M), H(N), HMAX(N), W(NW), X(N), Z(N)
      INTEGER           IG(NIG), IP(NIP), IRN(NIRN)
C     .. Scalars in Common ..
      INTEGER           L1, L2, N1, NC, NFUN, NG, NIMP
      LOGICAL           REPEAT
C     .. Local Scalars ..
      DOUBLE PRECISION  DER, HJL, HM, ONE, ROUND, TRUNC, XJ, ZERO
      INTEGER           I, ICC, J, K, K1, K2, L, MJ, N2, NC1
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AD02RA/NFUN, NIMP, REPEAT, L1, L2, NG, N1, NC
C     .. Data statements ..
      DATA              ZERO/0.D0/, ONE/1.D0/
C     .. Executable Statements ..
      IF (IND.GE.1 .AND. IND.LE.3) GO TO 20
      IND = 1
      GO TO 880
   20 GO TO (40,460,640) IND
   40 NFUN = 0
      N1 = N + 1
      IF (NIRN.GT.0 .AND. NIP.GE.N+1 .AND. NIG.GE.2*N+1 .AND. NW.GE.M+N)
     *    GO TO 60
      IND = 1
      GO TO 880
C     IG(N1+J),J=1,2,... ARE COLUMN NUMBERS FOR GROUP 1 THEN GROUP
C     2 ETC
   60 N2 = N - 1
      IF (N2.EQ.0) GO TO 100
      DO 80 I = 1, N2
         IF (IP(I+1).GE.IP(I)) GO TO 80
         IND = 2
         GO TO 880
   80 CONTINUE
  100 N2 = IP(N+1) - 1
      IF (N2.EQ.0) GO TO 140
      DO 120 I = 1, N2
         IF (IRN(I).GT.0 .AND. IRN(I).LE.M) GO TO 120
         IND = 2
         GO TO 880
  120 CONTINUE
  140 IF (IG(1).NE.0) GO TO 300
      ICC = 1
      DO 160 J = 1, N
         IG(J+1) = 0
         MJ = M + J
         W(MJ) = ZERO
  160 CONTINUE
C     W(M+J) IS ZERO IF COLUMN J HAS NOT BEEN INCLUDED IN A GROUP
C     YET
      DO 280 NC1 = 1, N1
         NC = NC1
         IG(NC) = ICC
         DO 180 I = 1, M
            W(I) = ZERO
  180    CONTINUE
C        W(I),I=1,M HOLDS THE BOOLEAN PATTERN OF THE UNION OF THE
C        COLUMNS SO FAR INCLUDED IN GROUP NC
         DO 260 J = 1, N
            MJ = M + J
            IF (W(MJ).NE.ZERO) GO TO 260
            K1 = IP(J)
            K2 = IP(J+1) - 1
            IF (K2.LT.K1) GO TO 260
            N2 = 0
            DO 200 K = K1, K2
               I = IRN(K)
               IF (I.GT.N2) GO TO 200
               IND = 2
               GO TO 880
  200       CONTINUE
            DO 220 K = K1, K2
               I = IRN(K)
               IF (W(I).NE.ZERO) GO TO 260
  220       CONTINUE
C           ACCEPT COLUMN
            DO 240 K = K1, K2
               I = IRN(K)
               W(I) = ONE
  240       CONTINUE
            I = N1 + ICC
            IG(I) = J
            ICC = ICC + 1
            W(MJ) = ONE
  260    CONTINUE
         IF (ICC.EQ.IG(NC)) GO TO 300
  280 CONTINUE
C     PRESERVE X IN W(M+1),...,W(M+N).   CHECK ALL H(J).
  300 DO 320 J = 1, N
         HM = -HMAX(1)
         IF (HM.LT.0.D0) HM = HMAX(J)
         XJ = ABS(X(J))
         IF (ADJUST) H(J) = (XJ+MAX(EPS*XJ,MIN(H(J),HM),EPS1*HM)) - XJ
         MJ = M + J
         W(MJ) = X(J)
  320 CONTINUE
C     Z(J) HOLDS MAXIMUM OF ESTIMATED RATIO OF TRUNCATION ERROR TO
C     ROUNDOFF ERROR IN COLUMN J,J=1,2,...,N. IP(J) IS NEGATED IF
C     H(J) HAS DECREASED DURING THIS CALL  OF TD02.
      NIMP = 0
  340 NIMP = NIMP + 1
      DO 360 J = 1, N
         Z(J) = 1.D0
  360 CONTINUE
C     FIND INITIAL APPROXIMATE DERIVATIVES
      NG = 0
  380 NG = NG + 1
      L1 = N1 + IG(NG)
      L2 = N1 + IG(NG+1) - 1
      IF (L2.LT.L1) GO TO 540
      DO 400 L = L1, L2
         J = IG(L)
         IF (H(J)) 400, 400, 420
  400 CONTINUE
      GO TO 520
  420 DO 440 L = L1, L2
         J = IG(L)
         H(J) = ABS(H(J))
         X(J) = X(J) + H(J)
  440 CONTINUE
      IND = 2
      RETURN
  460 NFUN = NFUN + 1
      DO 500 L = L1, L2
         J = IG(L)
         MJ = M + J
         X(J) = W(MJ)
         K1 = ABS(IP(J))
         K2 = ABS(IP(J+1)) - 1
         IF (K2.LT.K1) GO TO 500
         DO 480 K = K1, K2
            I = IRN(K)
            A(K) = (W(I)-F(I))/H(J)
  480    CONTINUE
  500 CONTINUE
  520 CONTINUE
      IF (NG.LT.N) GO TO 380
  540 IF ( .NOT. ADJUST) GO TO 860
C     ESTIMATE ERRORS AND IMPROVED STEP-LENGTHS
      NG = 0
  560 NG = NG + 1
      L1 = N1 + IG(NG)
      L2 = N1 + IG(NG+1) - 1
      IF (L2.LT.L1) GO TO 720
      DO 580 L = L1, L2
         J = IG(L)
         IF (H(J)) 580, 580, 600
  580 CONTINUE
      GO TO 700
  600 DO 620 L = L1, L2
         J = IG(L)
         X(J) = X(J) - H(J)
  620 CONTINUE
      IND = 3
      RETURN
  640 NFUN = NFUN + 1
      DO 680 L = L1, L2
         J = IG(L)
         MJ = M + J
         X(J) = W(MJ)
         K1 = ABS(IP(J))
         K2 = ABS(IP(J+1)) - 1
         IF (K2.LT.K1) GO TO 680
         DO 660 K = K1, K2
            I = IRN(K)
            DER = (F(I)-W(I))/H(J)
            ROUND = EPS*(0.5D0*(ABS(W(I))+ABS(A(K)*H(J)+F(I)))
     *              +MAX(ABS(A(K)),ABS(DER))*(ABS(X(J))+H(J)))/H(J)
            A(K) = (A(K)+DER)/2.D0
            TRUNC = DER - A(K)
            IF (ROUND.EQ.0.D0) ROUND = 1.D0
            Z(J) = MAX(Z(J),ABS(TRUNC/ROUND))
  660    CONTINUE
  680 CONTINUE
  700 CONTINUE
      IF (NG.LT.N) GO TO 560
C     FIND IMPROVED STEPS AND DECIDE WHETHER FURTHER SWEEP IS NEEDED
  720 REPEAT = .FALSE.
      DO 800 J = 1, N
         IF (H(J).LT.0.D0) GO TO 800
         HJL = H(J)
         XJ = ABS(X(J))
         HM = -HMAX(1)
         IF (HM.LT.0.D0) HM = HMAX(J)
         IF (Z(J).LT.UMIN) GO TO 760
         H(J) = (MAX(H(J)*SQRT(UAIM/Z(J)),EPS*XJ,EPS1*HM)+XJ) - XJ
         IF (H(J).LT.HJL) IP(J) = -IP(J)
         IF (Z(J).GT.UMAX) GO TO 780
  740    H(J) = -H(J)
         GO TO 800
  760    H(J) = (MIN(SQRT(UAIM/Z(J))*H(J),HM)+XJ) - XJ
  780    IF (ABS(HJL/H(J)-1.D0).LE.0.001D0) GO TO 740
         IF (NIMP.GE.6 .AND. IP(J).LT.0) GO TO 740
         REPEAT = .TRUE.
  800 CONTINUE
      IF ( .NOT. REPEAT) GO TO 820
      IF (NIMP.LT.60) GO TO 340
  820 DO 840 J = 1, N
         IP(J) = ABS(IP(J))
         H(J) = ABS(H(J))
  840 CONTINUE
  860 W(1) = NFUN
      W(2) = NC - 1
      IND = 0
  880 IFAIL = IND
      IND = 0
      RETURN
      END
