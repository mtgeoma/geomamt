      SUBROUTINE G01DAZ(N,NDIV2,CRLN,PROB1,PBEST,PP,ETOL,ERREST,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     G01DAZ COMPUTES NORMAL SCORES FOR G01DAF
C     PARAMETERS -
C     NDIV2     (N/2) + 1
C     CRLN,PROB1,PBEST  WORK ARRAYS OF DIMENSION NDIV2
C     N,PP,ETOL,ERREST,IFAIL  AS DESCRIBED IN G01DAF
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERREST, ETOL
      INTEGER           IFAIL, N, NDIV2
C     .. Array Arguments ..
      DOUBLE PRECISION  CRLN(NDIV2), PBEST(NDIV2), PP(N), PROB1(NDIV2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, A3, A4, A5, CK1, CK1A, CK2, CK2A, E,
     *                  EARG, ECONST, ELAST, EMAX, ENEG, EPS, EXP1LN,
     *                  EXP2LN, F1, F2, RIMIN1, RINC, RK, RLCK1, RLCK2,
     *                  RLMV1, RLMV2, RLV1, RLV2, RN, RNMINI, S, SMALL,
     *                  SRPLN, T3P, TWOPI, V1, V2, WMV2, X, X1, X2, XMIN
      INTEGER           I, IM1, J, K, KOUNT, LIMHI, LOWLIM, M, MFAIL
      LOGICAL           CANCOM
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02AKF, X02AMF
      EXTERNAL          X01AAF, X02AJF, X02AKF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, LOG, EXP, DBLE, SQRT
C     .. Executable Statements ..
      EPS = X02AJF()
      ENEG = LOG(X02AMF())
      XMIN = X02AKF()
      MFAIL = 0
      IF (N.GT.1) GO TO 40
      IF (N.EQ.1) GO TO 20
      MFAIL = 1
      GO TO 800
   20 PP(1) = 0.0D0
      ERREST = 0.0D0
      GO TO 800
   40 IF (ETOL.GT.0.0D0) GO TO 60
      MFAIL = 2
      GO TO 800
   60 SMALL = MAX(LOG(EPS)*1.5D0,ENEG)
      TWOPI = 2.0D0*X01AAF(TWOPI)
      SRPLN = LOG(TWOPI)*0.5D0
      M = N/2
      RN = N
      RK = RN
      CRLN(1) = LOG(RN) - SRPLN
      IF (M.LT.2) GO TO 100
      DO 80 I = 2, M
         RK = RK - 1.0D0
         IM1 = I - 1
         CRLN(I) = CRLN(IM1) + LOG(RK) - LOG(DBLE(IM1))
   80 CONTINUE
  100 RINC = 0.04D0
      IF (N.GE.300) RINC = 0.02D0
      IF (N.GE.1100) RINC = 0.01D0
      ELAST = ETOL + 2.0D0
      EMAX = ETOL + 1.0D0
      RINC = RINC/2.0D0
  120 IF (EMAX.LE.ETOL .OR. EMAX.GE.ELAST) GO TO 700
      ELAST = EMAX
      DO 140 J = 1, M
         PP(J) = 0.0D0
  140 CONTINUE
      T3P = RINC/(3.0D0*SQRT(TWOPI))
      X = 0.0D0
      S = 0.5D0
      A1 = 1.0D0
      F2 = 4.0D0*RINC/3.0D0
      F1 = F2 + F2
      ECONST = (RN-1.0D0)*LOG(0.5D0)
      DO 180 I = 1, M
         EARG = CRLN(I) + ECONST
         PROB1(I) = 0.0D0
         IF (EARG.LT.ENEG) GO TO 180
         IF (F2.LT.1.0D0) GO TO 160
         PROB1(I) = F2*EXP(EARG)
         GO TO 180
  160    IF (EXP(EARG).LT.XMIN/F2) GO TO 180
         PROB1(I) = F2*EXP(EARG)
  180 CONTINUE
      KOUNT = 0
      CANCOM = .TRUE.
      LIMHI = 1
      LOWLIM = M
      KOUNT = KOUNT + 1
  200 IF ( .NOT. (CANCOM .AND. (LIMHI.LE.M .OR. LOWLIM.GT.0)))
     *    GO TO 620
      X = X + RINC
      EARG = -0.5D0*X*X
      A2 = 0.0D0
      IF (EARG.GE.ENEG) A2 = EXP(EARG)
      X = X + RINC
      EARG = -0.5D0*X*X
      A3 = 0.0D0
      IF (EARG.GE.ENEG) A3 = EXP(EARG)
      X1 = X
      X = X + RINC
      EARG = -0.5D0*X*X
      A4 = 0.0D0
      IF (EARG.GE.ENEG) A4 = EXP(EARG)
      X = X + RINC
      EARG = -0.5D0*X*X
      A5 = 0.0D0
      IF (EARG.GE.ENEG) A5 = EXP(EARG)
      X2 = X
      V1 = S + (A1+4.0D0*A2+A3)*T3P
      V2 = V1 + (A3+4.0D0*A4+A5)*T3P
      S = V2
      A1 = A5
      WMV2 = 1.0D0 - V2
      IF (V2.LT.1.0D0 .AND. WMV2.GT.0.0D0) GO TO 220
      CANCOM = .FALSE.
      GO TO 600
  220 RLV1 = LOG(V1)
      RLV2 = LOG(V2)
      RLMV1 = LOG(1.0D0-V1)
      RLMV2 = LOG(WMV2)
      EXP1LN = -0.5D0*X1*X1
      EXP2LN = -0.5D0*X2*X2
      I = M
  240 IF (I.LT.LIMHI) GO TO 360
      RIMIN1 = I - 1
      RNMINI = N - I
      RLCK1 = EXP1LN + CRLN(I) + RIMIN1*RLV1 + RNMINI*RLMV1
      RLCK2 = EXP2LN + CRLN(I) + RIMIN1*RLV2 + RNMINI*RLMV2
      IF (RLCK1.GE.SMALL) GO TO 260
      LIMHI = I + 1
      GO TO 360
  260 CK1 = 0.0D0
      IF (RLCK1.LT.ENEG) GO TO 300
      IF (F1.LT.1.0D0) GO TO 280
      CK1 = EXP(RLCK1)*F1
      GO TO 300
  280 IF (EXP(RLCK1).LT.XMIN/F1) GO TO 300
      CK1 = EXP(RLCK1)*F1
  300 CK2 = 0.0D0
      IF (RLCK2.LT.ENEG) GO TO 340
      IF (F2.LT.1.0D0) GO TO 320
      CK2 = EXP(RLCK2)*F2
      GO TO 340
  320 IF (EXP(RLCK2).LT.XMIN/F2) GO TO 340
      CK2 = EXP(RLCK2)*F2
  340 PP(I) = PP(I) + CK1*X1 + CK2*X2
      PROB1(I) = PROB1(I) + CK1 + CK2
      I = I - 1
      GO TO 240
  360 I = 1
  380 IF (I.GT.LOWLIM) GO TO 600
      RIMIN1 = I - 1
      RNMINI = N - I
      RLCK1 = EXP1LN + CRLN(I) + RIMIN1*RLMV1 + RNMINI*RLV1
      RLCK2 = EXP2LN + CRLN(I) + RIMIN1*RLMV2 + RNMINI*RLV2
      IF (RLCK2.LT.RLCK1 .AND. RLCK1.LT.SMALL) GO TO 560
      CK1 = 0.0D0
      IF (RLCK1.LT.ENEG) GO TO 420
      IF (F1.LT.1.0D0) GO TO 400
      CK1 = EXP(RLCK1)*F1
      GO TO 420
  400 IF (EXP(RLCK1).LT.XMIN/F1) GO TO 420
      CK1 = EXP(RLCK1)*F1
  420 CK2 = 0.0D0
      IF (RLCK2.LT.ENEG) GO TO 460
      IF (F2.LT.1.0D0) GO TO 440
      CK2 = EXP(RLCK2)*F2
      GO TO 460
  440 IF (EXP(RLCK2).LT.XMIN/F2) GO TO 460
      CK2 = EXP(RLCK2)*F2
  460 CK1A = 0.0D0
      IF (X1.LT.1.0D0) GO TO 480
      CK1A = CK1*X1
      GO TO 500
  480 IF (CK1.GE.XMIN/X1) CK1A = CK1*X1
  500 CK2A = 0.0D0
      IF (X2.LT.1.0D0) GO TO 520
      CK2A = CK2*X2
      GO TO 540
  520 IF (CK2.GE.XMIN/X2) CK2A = CK2*X2
  540 PP(I) = PP(I) - CK1A - CK2A
      PROB1(I) = PROB1(I) + CK1 + CK2
      GO TO 580
  560 LOWLIM = RIMIN1
  580 I = I + 1
      GO TO 380
  600 KOUNT = KOUNT + 1
      GO TO 200
  620 EMAX = 0.0D0
      DO 640 I = 1, M
         E = ABS(1.0D0-PROB1(I))
         IF (E.GT.EMAX) EMAX = E
  640 CONTINUE
      EMAX = EMAX*X
      IF (EMAX.LE.ETOL) GO TO 680
      DO 660 I = 1, M
         PBEST(I) = PP(I)
  660 CONTINUE
  680 RINC = RINC/2.0D0
      GO TO 120
  700 IF (EMAX.LE.ELAST) GO TO 740
      EMAX = ELAST
      DO 720 I = 1, M
         PP(I) = PBEST(I)
  720 CONTINUE
  740 I = M + 1
      IF (M+M.EQ.N) GO TO 760
      PP(M+1) = 0.0D0
      I = I + 1
  760 J = M
      DO 780 K = I, N
         PP(K) = -PP(J)
         J = J - 1
  780 CONTINUE
      ERREST = EMAX
      IF (EMAX.GT.ETOL) MFAIL = 3
  800 IF (MFAIL.NE.0) GO TO 820
      IFAIL = 0
      RETURN
  820 IFAIL = MFAIL
      RETURN
      END
