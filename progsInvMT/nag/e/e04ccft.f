      SUBROUTINE E04CCF(N,X,FMIN,EPS,N1,PDSTAR,PSTAR,PBAR,STEP,Y,P,
     *                  FUNCT,MONIT,MAXIT,IFAIL)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED.
C     MARK 4.5 REISSUE. LER-F7
C     MARK 8 REVISED. IER-220 (MAR 1980)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-716 (DEC 1989).
C     MARK 16A REVISED. IER-984 (JUN 1993).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04CCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, FMIN
      INTEGER           IFAIL, MAXIT, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  P(N1,N), PBAR(N), PDSTAR(N), PSTAR(N), STEP(N),
     *                  X(N), Y(N1)
C     .. Subroutine Arguments ..
      EXTERNAL          FUNCT, MONIT
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C, CENT, COEFF, DERIV, DERIV2, F1, F2, F3,
     *                  FMAX, R, SERROR, X1, X2, X3, XMIN, YDSTAR,
     *                  YMEAN, YSTAR, F1FIX
      INTEGER           H, I, IV, J, K, L, LASTMX, MCOUNT, NCALL, NP1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT, DBLE
C     .. Executable Statements ..
      FMAX = 1000.D0
      FMIN = 0.D0
      NCALL = 0
      R = 0
      IF (N.LT.1 .OR. EPS.LT.X02AJF()
     *    .OR. N1.NE.N+1 .OR. MAXIT.LT.1) GO TO 1220
      CALL FUNCT(N,X,F1FIX)
      NCALL = NCALL + 1
   20 DO 300 I = 1, N
         F1 = F1FIX
         X1 = 0.D0
         COEFF = 1.D0
         DO 80 J = 1, N
            IF (I-J) 60, 40, 60
   40       PSTAR(J) = X(J) + COEFF
            GO TO 80
   60       PSTAR(J) = X(J)
   80    CONTINUE
         CALL FUNCT(N,PSTAR,F2)
         NCALL = NCALL + 1
         X2 = 1.D0
         PSTAR(I) = PSTAR(I) + COEFF
         CALL FUNCT(N,PSTAR,F3)
         NCALL = NCALL + 1
         X3 = 2.D0
  100    IF (NCALL.GT.MAXIT) GO TO 1120
         DERIV = (X2-X3)*F1 + (X3-X1)*F2 + (X1-X2)*F3
         IF (ABS(DERIV)-X02AJF()) 140, 120, 120
  120    DERIV2 = DERIV/(X1-X2)/(X2-X3)/(X3-X1)
         IF (DERIV2) 240, 140, 140
  140    IF (F1-F3) 160, 200, 200
  160    IF (X1.LE.-5.0D0) GO TO 180
         F3 = F2
         X3 = X2
         F2 = F1
         X2 = X1
         X1 = X1 - COEFF
         PSTAR(I) = X(I) + X1
         CALL FUNCT(N,PSTAR,F1)
         NCALL = NCALL + 1
         GO TO 100
  180    XMIN = -5.0D0
         GO TO 280
  200    IF (X3.GE.5.0D0) GO TO 220
         F1 = F2
         X1 = X2
         F2 = F3
         X2 = X3
         X3 = X3 + COEFF
         PSTAR(I) = X(I) + X3
         CALL FUNCT(N,PSTAR,F3)
         NCALL = NCALL + 1
         GO TO 100
  220    XMIN = 5.0D0
         GO TO 280
  240    XMIN = .5D0*((X2**2-X3**2)*F1+(X3**2-X1**2)*F2+(X1**2-X2**2)
     *          *F3)/DERIV
         IF (XMIN.NE.0.0D0) GO TO 260
         XMIN = 0.1D0
         GO TO 280
  260    IF (ABS(XMIN).LT.0.1D0) XMIN = SIGN(0.1D0,XMIN)
         IF (ABS(XMIN).GT.5.0D0) XMIN = SIGN(5.0D0,XMIN)
  280    STEP(I) = XMIN
  300 CONTINUE
      NP1 = N + 1
      DO 380 I = 1, NP1
         DO 360 J = 1, N
            IF (I-J-1) 320, 340, 320
  320       PSTAR(J) = X(J)
            P(I,J) = X(J)
            GO TO 360
  340       PSTAR(J) = X(J) + STEP(J)
            P(I,J) = X(J) + STEP(J)
  360    CONTINUE
         NCALL = NCALL + 1
         CALL FUNCT(N,PSTAR,Y(I))
  380 CONTINUE
      A = 1.D0
      B = .5D0
      C = 2.D0
      LASTMX = 0
      MCOUNT = 0
      K = 0
  400 K = K + 1
      FMAX = Y(1)
      FMIN = Y(1)
      H = 1
      L = 1
      DO 480 I = 2, N1
         IF (Y(I)-FMAX) 440, 440, 420
  420    FMAX = Y(I)
         H = I
         GO TO 480
  440    IF (Y(I)-FMIN) 460, 480, 480
  460    FMIN = Y(I)
         L = I
  480 CONTINUE
      IF (LASTMX-H) 580, 500, 580
  500 MCOUNT = MCOUNT + 1
      IF (MCOUNT-5) 600, 520, 600
  520 IF (H-1) 560, 540, 560
  540 H = 2
      FMAX = Y(H)
      GO TO 600
  560 H = 1
      FMAX = Y(H)
      GO TO 600
  580 LASTMX = H
      MCOUNT = 0
  600 CONTINUE
      CALL MONIT(FMIN,FMAX,P,N,N1,NCALL)
      IF (K.EQ.1) GO TO 620
      IF (SERROR-EPS) 1140, 620, 620
  620 IF (NCALL-MAXIT) 640, 1160, 1160
  640 DO 700 J = 1, N
         CENT = 0.D0
         DO 680 I = 1, N1
            IF (I-H) 660, 680, 660
  660       CENT = CENT + P(I,J)
  680    CONTINUE
         PBAR(J) = CENT/DBLE(N)
  700 CONTINUE
C     REFLECTION
      DO 720 I = 1, N
         PSTAR(I) = (1.D0+A)*PBAR(I) - A*P(H,I)
  720 CONTINUE
      CALL FUNCT(N,PSTAR,YSTAR)
      NCALL = NCALL + 1
      IF (YSTAR-FMIN) 740, 780, 780
C     EXPANSION
  740 DO 760 I = 1, N
         PDSTAR(I) = (1.D0+C)*PSTAR(I) - C*PBAR(I)
  760 CONTINUE
      CALL FUNCT(N,PDSTAR,YDSTAR)
      NCALL = NCALL + 1
      IF (YDSTAR-YSTAR) 980, 1020, 1020
C     CONTRACTION
  780 DO 820 I = 1, N
         IF (I-H) 800, 820, 800
  800    IF (YSTAR-Y(I)) 1020, 820, 820
  820 CONTINUE
      IF (FMAX-YSTAR) 880, 840, 840
  840 DO 860 I = 1, N
         P(H,I) = PSTAR(I)
  860 CONTINUE
  880 CONTINUE
      DO 900 I = 1, N
         PDSTAR(I) = B*P(H,I) + (1.D0-B)*PBAR(I)
  900 CONTINUE
      CALL FUNCT(N,PDSTAR,YDSTAR)
      NCALL = NCALL + 1
      IF (YDSTAR-FMAX) 980, 980, 920
  920 DO 960 I = 1, N1
         DO 940 J = 1, N
            PBAR(J) = (P(I,J)+P(L,J))*0.5D0
            P(I,J) = PBAR(J)
  940    CONTINUE
         CALL FUNCT(N,PBAR,Y(I))
         NCALL = NCALL + 1
  960 CONTINUE
      GO TO 1060
  980 DO 1000 J = 1, N
         P(H,J) = PDSTAR(J)
 1000 CONTINUE
      Y(H) = YDSTAR
      GO TO 1060
 1020 DO 1040 J = 1, N
         P(H,J) = PSTAR(J)
 1040 CONTINUE
      Y(H) = YSTAR
 1060 YMEAN = 0.D0
      SERROR = 0.D0
      DO 1080 I = 1, N1
         YMEAN = YMEAN + Y(I)
 1080 CONTINUE
      YMEAN = YMEAN/DBLE(N+1)
      DO 1100 I = 1, N1
         SERROR = SERROR + (Y(I)-YMEAN)**2
 1100 CONTINUE
      SERROR = SQRT(SERROR/DBLE(N+1))
      GO TO 400
 1120 IV = 2
      GO TO 1240
 1140 IV = 0
      GO TO 1180
 1160 IV = 2
 1180 DO 1200 I = 1, N
         X(I) = P(L,I)
 1200 CONTINUE
      GO TO 1240
 1220 IV = 1
 1240 IFAIL = P01ABF(IFAIL,IV,SRNAME,0,P01REC)
      RETURN
      END
