      SUBROUTINE D02AGF(H,ERROR,PARERR,PARAM,C,N,N1,M1,AUX,BCAUX,RAAUX,
     *                  PRSOL,MAT,COPY,WSPACE,WSPAC1,WSPAC2,IFAIL)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 15 REVISED. IER-898 (APR 1991). (LAPACK)
C
C     SOLVES A GENERAL BOUNDARY VALUE
C     PROBLEM FOR N DIFFERENTIAL EQUATIONS
C     IN N1 PARAMETERS USING A SHOOTING
C     AND MATCHING TECHNIQUE.  EPS IS THE
C     LARGEST REAL VARIABLE SUCH THAT 1+EPS=1
C     ALL IMPLICITLY DECLARED REALS MAY BE USED DOUBLE-LENGTH
C     THE ARRAY COPY IS REDUNDANT
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02AGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H
      INTEGER           IFAIL, M1, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  C(M1,N), COPY(N1,N1), ERROR(N), MAT(N1,N1),
     *                  PARAM(N1), PARERR(N1), WSPAC1(N), WSPAC2(N),
     *                  WSPACE(N,9)
C     .. Subroutine Arguments ..
      EXTERNAL          AUX, BCAUX, PRSOL, RAAUX
C     .. Local Scalars ..
      DOUBLE PRECISION  C1, D, DIST, DUM, EPS, H0, OLDRES, PERT, R,
     *                  RESID, X, X1
      INTEGER           COUNT, COUNT1, CT, EM, I, ITEST, J, K, M, ONE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02AGZ, F07ADG, F07AEG
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      EPS = X02AJF()
      M = M1 - 1
      IF (N1.LE.N) GO TO 20
      EM = 1
      GO TO 940
   20 COUNT = 0
      H0 = H
      ONE = 1
      EM = -1
C
C     FORMS THE RESIDUALS AT THE
C     MATCHING POINT
   40 CALL RAAUX(X,X1,R,PARAM)
      IF ((X-R)*(X1-R).LE.0.0D0) GO TO 60
      EM = 3
      GO TO 940
   60 IF (H0*(X1-X).LT.0.0D0) H0 = -H0
      CALL BCAUX(WSPAC1,WSPAC2,PARAM)
      H = H0
      DO 80 I = 1, N
         WSPACE(I,8) = WSPAC2(I)
   80 CONTINUE
      I = 1
      CALL D02AGZ(X,WSPAC1,ERROR,ONE,N,N1,I,R-X,H,AUX,WSPACE,WSPAC2,
     *            PARAM)
      IF (I.EQ.0) GO TO 100
      EM = 4
      GO TO 940
  100 DO 120 I = 1, N
         DUM = WSPACE(I,8)
         WSPACE(I,8) = -WSPAC1(I)
         WSPAC1(I) = DUM
  120 CONTINUE
      H = -H0
      I = 1
      CALL D02AGZ(X1,WSPAC1,ERROR,ONE,N,N1,I,R-X1,H,AUX,WSPACE,WSPAC2,
     *            PARAM)
      IF (I.EQ.0) GO TO 140
      EM = 4
      GO TO 940
  140 RESID = 0.0D0
      CT = 0
      DO 160 I = 1, N1
         D = WSPAC1(I)
         DUM = WSPACE(I,8)
         WSPACE(I,8) = D + DUM
         D = 1.0D0 + ABS(DUM) + ABS(D)
         DUM = WSPACE(I,8)
         IF (ABS(DUM).LT.ERROR(I)*D) CT = CT + 1
         WSPAC1(I) = DUM
         RESID = RESID + DUM*DUM
  160 CONTINUE
      CALL PRSOL(PARAM,RESID,N1,WSPAC1)
      IF (EM.NE.-1) GO TO 620
  180 COUNT = COUNT + 1
      IF (COUNT.NE.12) GO TO 200
      EM = 7
      GO TO 940
C
C     FORMS THE JACOBIAN BY NUMERICAL
C     DIFFERENTIATION
  200 DO 360 K = 1, N1
         PERT = 10.0D0*PARERR(K)*(1.0D0+ABS(PARAM(K)))
         PARAM(K) = PERT + PARAM(K)
         CALL RAAUX(X,X1,R,PARAM)
         IF ((X-R)*(X1-R).LE.0.0D0) GO TO 220
         EM = 3
         GO TO 940
  220    IF (H0*(X1-X).LT.0.0D0) H0 = -H0
         H = H0
         CALL BCAUX(WSPAC1,WSPAC2,PARAM)
         DO 240 I = 1, N
            WSPACE(I,7) = WSPAC2(I)
  240    CONTINUE
         I = 1
         CALL D02AGZ(X,WSPAC1,ERROR,ONE,N,N1,I,R-X,H,AUX,WSPACE,WSPAC2,
     *               PARAM)
         IF (I.EQ.0) GO TO 260
         EM = 2
         GO TO 940
  260    DO 280 I = 1, N1
            MAT(I,K) = WSPAC1(I)
  280    CONTINUE
         H = -H0
         I = 1
         DO 300 I = 1, N
            WSPAC1(I) = WSPACE(I,7)
  300    CONTINUE
         CALL D02AGZ(X1,WSPAC1,ERROR,ONE,N,N1,I,R-X1,H,AUX,WSPACE,
     *               WSPAC2,PARAM)
         IF (I.EQ.0) GO TO 320
         EM = 2
         GO TO 940
  320    DO 340 I = 1, N1
            MAT(I,K) = (MAT(I,K)-WSPAC1(I)+WSPACE(I,8))/PERT
            IF (ABS(MAT(I,K)).LT.5.0D0*EPS*ABS(WSPACE(I,8))/PERT)
     *          MAT(I,K) = 0.0D0
  340    CONTINUE
         PARAM(K) = PARAM(K) - PERT
  360 CONTINUE
      ITEST = 1
      EM = -3
C
C     PERFORMS COLUMN SCALING ON THE JACOBIAN
C     AND FORMS A TRIANGULAR DECOMPOSITION
      DO 440 I = 1, N1
         D = 0.0D0
         DO 380 J = 1, N1
            IF (ABS(MAT(J,I)).GT.D) D = ABS(MAT(J,I))
  380    CONTINUE
         IF (D.NE.0.0D0) GO TO 400
         EM = 5
         GO TO 940
  400    DO 420 J = 1, N1
            MAT(J,I) = MAT(J,I)/D
  420    CONTINUE
         WSPACE(I,7) = D
  440 CONTINUE
      CALL F07ADG(N1,N1,MAT,N1,WSPAC1,I)
      IF (I.EQ.0) GO TO 460
      EM = 5
      GO TO 940
  460 DO 480 I = 1, N1
         WSPACE(I,6) = WSPAC1(I)
  480 CONTINUE
C
C     USES A GENERALISED NEWTON RAPHSON
C     TECHNIQUE TO SOLVE THE NONLINEAR
C     EQUATIONS AT THE MATCHING POINT
  500 OLDRES = RESID
      COUNT1 = 0
      DO 520 I = 1, N1
         WSPAC1(I) = WSPACE(I,6)
         WSPACE(I,1) = WSPACE(I,8)
  520 CONTINUE
      CALL F07AEG('No transpose',N1,ONE,MAT,N1,WSPAC1,WSPACE,N,I)
      DO 540 I = 1, N1
         WSPACE(I,1) = WSPACE(I,1)/WSPACE(I,7)
         PARAM(I) = PARAM(I) + WSPACE(I,1)
  540 CONTINUE
      IF (CT.LT.N1) GO TO 580
      DO 560 I = 1, N1
         IF (ABS(WSPACE(I,1)).GT.PARERR(I)*(1.0D0+ABS(PARAM(I))))
     *       GO TO 580
  560 CONTINUE
      EM = -5
      GO TO 740
  580 DO 600 I = 1, N1
         WSPACE(I,1) = -WSPACE(I,1)
  600 CONTINUE
      GO TO 40
  620 IF (COUNT1.NE.0) GO TO 640
      IF (RESID.GE.OLDRES/10.0D0) GO TO 640
      EM = -2
      ITEST = 0
      GO TO 500
  640 IF (RESID.LT.OLDRES) GO TO 180
      IF (COUNT1.NE.3) GO TO 700
      IF (ITEST.EQ.0) GO TO 660
      EM = 6
      GO TO 940
  660 DO 680 I = 1, N1
         PARAM(I) = PARAM(I) + WSPACE(I,1)
  680 CONTINUE
      EM = -1
      GO TO 40
  700 COUNT1 = COUNT1 + 1
      EM = -4
      DO 720 I = 1, N1
         WSPACE(I,1) = WSPACE(I,1)/2.0D0
         PARAM(I) = PARAM(I) + WSPACE(I,1)
  720 CONTINUE
      GO TO 40
C
C     CALCULATES THE FINAL SOLUTION
  740 IF (M.LE.0) GO TO 940
      CALL RAAUX(X,X1,R,PARAM)
      IF ((X-R)*(X1-R).LE.0.0D0) GO TO 760
      EM = 3
      GO TO 940
  760 IF (H0*(X1-X).LT.0.0D0) H0 = -H0
      H = H0
      CALL BCAUX(WSPAC1,WSPAC2,PARAM)
      DO 780 I = 1, N
         WSPACE(I,7) = WSPAC2(I)
  780 CONTINUE
      DIST = (X1-X)/DBLE(M)
      J = 1
      C1 = X
      K = 1
  800 DO 820 I = 1, N
         C(J,I) = WSPAC1(I)
  820 CONTINUE
      IF ((R-C1-0.25D0*DIST)*DIST.LE.0.0D0) GO TO 860
      I = 1
      CALL D02AGZ(C1,WSPAC1,ERROR,ONE,N,N1,I,DIST,H,AUX,WSPACE,WSPAC2,
     *            PARAM)
      IF (I.EQ.0) GO TO 840
      EM = 4
      GO TO 940
  840 J = J + K
      GO TO 800
  860 IF (K.EQ.-1) GO TO 900
      DIST = -DIST
      C1 = X1
      H = -H0
      J = M1
      K = -1
      DO 880 I = 1, N
         WSPAC1(I) = WSPACE(I,7)
  880 CONTINUE
      GO TO 800
  900 CALL BCAUX(WSPAC1,WSPAC2,PARAM)
      DO 920 I = 1, N
         C(1,I) = WSPAC1(I)
         C(M1,I) = WSPAC2(I)
  920 CONTINUE
  940 IF (EM.LE.0) IFAIL = 0
      IF (EM.GT.0) IFAIL = P01ABF(IFAIL,EM,SRNAME,0,P01REC)
      RETURN
C     END OF D02AGF
      END
