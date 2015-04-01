      SUBROUTINE F02AKZ(N,LOW,IUPP,ACHEPS,HR,IHR,HI,IHI,WR,WI,VR,IVR,VI,
     *                  IVI,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACHEPS
      INTEGER           IERR, IHI, IHR, IUPP, IVI, IVR, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  HI(IHI,N), HR(IHR,N), VI(IVI,N), VR(IVR,N),
     *                  WI(N), WR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AAHR, AHR, ANORM, SI, SR, TI, TR, XI, XR, XXI,
     *                  XXR, YI, YR, YYI, YYR, ZI, ZR
      INTEGER           I, I1, IEN, IEN1, IIEN, ITN, ITS, IUPP1, IUPP11,
     *                  J, J1, JJ, K, K1, KK, LOW1, LOW11, M, M1, MM
C     .. External Subroutines ..
      EXTERNAL          A02AAF, A02ACF, F01AMZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
      LOW1 = LOW + 1
      IUPP1 = IUPP - 1
      LOW11 = LOW - 1
      DO 20 I = 1, LOW11
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
   20 CONTINUE
      IUPP11 = IUPP + 1
      DO 40 I = IUPP11, N
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
   40 CONTINUE
C     END OF ISOLATED ROOTS
      IEN = IUPP
   60 IF (IEN.LT.LOW) GO TO 600
      ITS = 0
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
   80 IF (LOW1.GT.IEN) GO TO 120
      DO 100 KK = LOW1, IEN
         K = LOW1 + IEN - KK
         K1 = K - 1
         AHR = ABS(HR(K,K1)) + ABS(HI(K,K1))
         AAHR = ACHEPS*(ABS(HR(K1,K1))+ABS(HI(K1,K1))+ABS(HR(K,K))
     *          +ABS(HI(K,K)))
         IF (AHR.LE.AAHR) GO TO 140
  100 CONTINUE
  120 K = LOW
  140 IF (K.EQ.IEN) GO TO 580
      IF (ITN.LE.0) GO TO 860
C     FORM SHIFT
      IF (ITS.EQ.10 .OR. ITS.EQ.20) GO TO 180
      SR = HR(IEN,IEN)
      SI = HI(IEN,IEN)
      IEN1 = IEN - 1
      XR = HR(IEN1,IEN)*HR(IEN,IEN1) - HI(IEN1,IEN)*HI(IEN,IEN1)
      XI = HR(IEN1,IEN)*HI(IEN,IEN1) + HI(IEN1,IEN)*HR(IEN,IEN1)
      IF (XR.EQ.0.0D0 .AND. XI.EQ.0.0D0) GO TO 200
      YR = (HR(IEN1,IEN1)-SR)/2.0D0
      YI = (HI(IEN1,IEN1)-SI)/2.0D0
      CALL A02AAF(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZR,ZI)
      IF ((YR*ZR+YI*ZI).GE.0.0D0) GO TO 160
      ZR = -ZR
      ZI = -ZI
  160 YYI = YI + ZI
      YYR = YR + ZR
      XXR = XR
      XXI = XI
      CALL A02ACF(XXR,XXI,YYR,YYI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 200
  180 IEN1 = IEN - 1
      SR = ABS(HR(IEN,IEN1)) + ABS(HR(IEN1,IEN-2))
      SI = ABS(HI(IEN,IEN1)) + ABS(HI(IEN1,IEN-2))
  200 IF (LOW.GT.IEN) GO TO 240
      DO 220 I = LOW, IEN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  220 CONTINUE
  240 TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
      J = K + 1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
      XR = ABS(HR(IEN1,IEN1)) + ABS(HI(IEN1,IEN1))
      YR = ABS(HR(IEN,IEN1)) + ABS(HI(IEN,IEN1))
      ZR = ABS(HR(IEN,IEN)) + ABS(HI(IEN,IEN))
      IF (J.GT.IEN1) GO TO 280
      DO 260 MM = J, IEN1
         M = J + IEN1 - MM
         YI = YR
         M1 = M - 1
         YR = ABS(HR(M,M1)) + ABS(HI(M,M1))
         XI = ZR
         ZR = XR
         XR = ABS(HR(M1,M1)) + ABS(HI(M1,M1))
         IF (YR.LE.(ACHEPS*ZR/YI*(ZR+XR+XI))) GO TO 300
  260 CONTINUE
  280 M = K
C     TRIANGULAR DECOMPOSITION H = L*R
  300 M1 = M + 1
      IF (M1.GT.IEN) GO TO 440
      DO 420 I = M1, IEN
         I1 = I - 1
         XR = HR(I1,I1)
         XI = HI(I1,I1)
         YR = HR(I,I1)
         YI = HI(I,I1)
         IF ((ABS(XR)+ABS(XI)).LT.(ABS(YR)+ABS(YI))) GO TO 320
         CALL A02ACF(YR,YI,XR,XI,ZR,ZI)
         WR(I) = -1.0D0
         GO TO 380
  320    IF (I1.GT.N) GO TO 360
C        INTERCHANGE ROWS OF HR AND HI
         DO 340 J = I1, N
            ZR = HR(I1,J)
            HR(I1,J) = HR(I,J)
            HR(I,J) = ZR
            ZI = HI(I1,J)
            HI(I1,J) = HI(I,J)
            HI(I,J) = ZI
  340    CONTINUE
  360    CALL A02ACF(XR,XI,YR,YI,ZR,ZI)
         WR(I) = 1.0D0
  380    HR(I,I1) = ZR
         HI(I,I1) = ZI
         IF (I.GT.N) GO TO 420
         DO 400 J = I, N
            HR(I,J) = HR(I,J) - ZR*HR(I1,J) + ZI*HI(I1,J)
            HI(I,J) = HI(I,J) - ZR*HI(I1,J) - ZI*HR(I1,J)
  400    CONTINUE
  420 CONTINUE
  440 IF (M1.GT.IEN) GO TO 80
C     COMPOSITION R*L = H
      DO 560 J = M1, IEN
         J1 = J - 1
         XR = HR(J,J1)
         XI = HI(J,J1)
         HR(J,J1) = 0.0D0
         HI(J,J1) = 0.0D0
C        INTERCHANGE COLUMNS OF HR, HI, VR, AND VI, IF NECESSARY
         IF (WR(J).LE.0.0D0) GO TO 500
         DO 460 I = 1, J
            ZR = HR(I,J1)
            HR(I,J1) = HR(I,J)
            HR(I,J) = ZR
            ZI = HI(I,J1)
            HI(I,J1) = HI(I,J)
            HI(I,J) = ZI
  460    CONTINUE
         IF (LOW.GT.IUPP) GO TO 500
         DO 480 I = LOW, IUPP
            ZR = VR(I,J1)
            VR(I,J1) = VR(I,J)
            VR(I,J) = ZR
            ZI = VI(I,J1)
            VI(I,J1) = VI(I,J)
            VI(I,J) = ZI
  480    CONTINUE
C        END INTERCHANGE COLUMNS
  500    DO 520 I = 1, J
            HR(I,J1) = HR(I,J1) + XR*HR(I,J) - XI*HI(I,J)
            HI(I,J1) = HI(I,J1) + XR*HI(I,J) + XI*HR(I,J)
  520    CONTINUE
         IF (LOW.GT.IUPP) GO TO 560
         DO 540 I = LOW, IUPP
            VR(I,J1) = VR(I,J1) + XR*VR(I,J) - XI*VI(I,J)
            VI(I,J1) = VI(I,J1) + XR*VI(I,J) + XI*VR(I,J)
  540    CONTINUE
C        END ACCUMULATE TRANSFORMATIONS
  560 CONTINUE
      GO TO 80
C     A ROOT FOUND
  580 WR(IEN) = HR(IEN,IEN) + TR
      WI(IEN) = HI(IEN,IEN) + TI
      IEN = IEN - 1
      GO TO 60
C     ALL ROOTS FOUND
  600 ANORM = 0.0D0
      DO 640 I = 1, N
         ANORM = ANORM + ABS(WR(I)) + ABS(WI(I))
         I1 = I + 1
         IF (I1.GT.N) GO TO 640
         DO 620 J = I1, N
            ANORM = ANORM + ABS(HR(I,J)) + ABS(HI(I,J))
  620    CONTINUE
  640 CONTINUE
      IF (ANORM.EQ.0.0D0) GO TO 840
C     BACKSUBSTITUTE TO FIND VECTORS OF UPPER TRIANGULAR FORM
      IF (N.LT.2) GO TO 720
      DO 700 IIEN = 2, N
         IEN = 2 + N - IIEN
         XR = WR(IEN)
         XI = WI(IEN)
         IEN1 = IEN - 1
         DO 680 JJ = 1, IEN1
            J = 1 + IEN1 - JJ
            YR = XR - WR(J)
            YI = XI - WI(J)
            IF (YR.EQ.0.0D0 .AND. YI.EQ.0.0D0) YR = ACHEPS*ANORM
            CALL A02ACF(HR(J,IEN),HI(J,IEN),YR,YI,ZR,ZI)
            HR(J,IEN) = ZR
            HI(J,IEN) = ZI
            J1 = J - 1
            IF (J1.EQ.0) GO TO 680
            DO 660 I = 1, J1
               HR(I,IEN) = HR(I,IEN) + HR(I,J)*ZR - HI(I,J)*ZI
               HI(I,IEN) = HI(I,IEN) + HR(I,J)*ZI + HI(I,J)*ZR
  660       CONTINUE
  680    CONTINUE
  700 CONTINUE
C     END BACKSUB
C     MULTIPLY BY TRANSFORMATION MATRIX TO GIVE VECTORS OF
C     ORIGINAL FULL MATRIX
  720 IUPP1 = IUPP + 1
      DO 800 J = 1, N
         M = MIN(J,LOW) - 1
         IF (M.LT.1) GO TO 760
         DO 740 I = 1, M
            VR(I,J) = HR(I,J)
            VI(I,J) = HI(I,J)
  740    CONTINUE
  760    M = J - 1
         IF (IUPP1.GT.M) GO TO 800
         DO 780 I = IUPP1, M
            VR(I,J) = HR(I,J)
            VI(I,J) = HI(I,J)
  780    CONTINUE
  800 CONTINUE
C     END VECTORS OF ISOLATED ROOTS
      IF (LOW.GT.N) GO TO 840
      DO 820 JJ = LOW, N
         J = LOW + N - JJ
         M = MIN(J-1,IUPP)
         IF (LOW.LE.M) CALL F01AMZ(VR(LOW,LOW),IVR,VI(LOW,LOW),IVI,
     *                             IUPP-LOW+1,M-LOW+1,HR(LOW,J),1,
     *                             HI(LOW,J),1,VR(LOW,J),VI(LOW,J))
  820 CONTINUE
  840 IERR = 0
      RETURN
  860 IERR = 1
      RETURN
      END
