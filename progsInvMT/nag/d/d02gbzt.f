      SUBROUTINE D02GBZ(E,F,N,C,D,IND,W)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     CHECKS FOR RANK-DEFICIENT BOUNDARY CONDITIONS
C     .. Scalar Arguments ..
      INTEGER           IND, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N,N), D(N,N), E(N,N), F(N,N), W(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, SMAX, X
      INTEGER           I, I1, J, K, M
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, DBLE
C     .. Executable Statements ..
      IND = 0
      EPS = 10.D0*X02AJF()
      DO 60 I = 1, N
         W(I) = DBLE(I)
         SMAX = 0.D0
         DO 20 J = 1, N
            SMAX = MAX(ABS(SMAX),ABS(E(I,J)),ABS(F(I,J)))
   20    CONTINUE
         DO 40 J = 1, N
            C(I,J) = E(I,J)/SMAX
            D(I,J) = F(I,J)/SMAX
   40    CONTINUE
   60 CONTINUE
      DO 180 I = 1, N
         SMAX = 0.D0
         DO 80 J = 1, N
            IF (ABS(C(J,I)).LE.SMAX .OR. W(J).LT.0.D0) GO TO 80
            M = J
            SMAX = ABS(C(J,I))
   80    CONTINUE
         IF (SMAX.EQ.0.D0) GO TO 180
         W(M) = -W(M)
         DO 160 J = 1, N
            IF (C(J,I).EQ.0.D0 .OR. W(J).LT.0.D0) GO TO 160
            X = C(J,I)/SMAX
            C(J,I) = 0.D0
            IF (I.EQ.N) GO TO 120
            I1 = I + 1
            DO 100 K = I1, N
               C(J,K) = C(J,K) - X*C(M,K)
  100       CONTINUE
  120       DO 140 K = 1, N
               D(J,K) = D(J,K) - X*D(M,K)
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
      DO 200 I = 1, N
         IF (W(I).GT.0.D0) GO TO 220
  200 CONTINUE
      GO TO 360
  220 DO 300 I = 1, N
         SMAX = 0.D0
         DO 240 J = 1, N
            IF (ABS(D(J,I)).LE.SMAX .OR. W(J).LT.0.D0) GO TO 240
            M = J
            SMAX = ABS(D(J,I))
  240    CONTINUE
         IF (SMAX.EQ.0.D0) GO TO 300
         W(M) = -W(M)
         DO 280 J = 1, N
            IF (D(J,I).EQ.0.D0 .OR. W(J).LT.0.D0) GO TO 280
            X = D(J,I)/SMAX
            D(J,I) = 0.D0
            IF (I.EQ.N) GO TO 280
            I1 = I + 1
            DO 260 K = I1, N
               D(J,K) = D(J,K) - X*D(M,K)
  260       CONTINUE
  280    CONTINUE
  300 CONTINUE
      DO 320 I = 1, N
         IF (W(I).GT.0.D0) GO TO 340
  320 CONTINUE
      GO TO 360
  340 IND = I
      RETURN
  360 DO 400 I = 1, N
         DO 380 J = 1, N
            IF (ABS(C(I,J)).GE.EPS .OR. ABS(D(I,J)).GE.EPS) GO TO 400
  380    CONTINUE
         GO TO 420
  400 CONTINUE
      RETURN
  420 IND = -I
      RETURN
      END
