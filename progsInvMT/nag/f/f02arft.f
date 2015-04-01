      SUBROUTINE F02ARF(N,LOW,IUPP,ACHEPS,INTGER,HR,IHR,HI,IHI,WR,WI,VR,
     *                  IVR,VI,IVI,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 6 REVISED  IER-83
C     MARK 9 REVISED. IER-326 (SEP 1981).
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMLR2
C     FINDS THE EIGENVALUES AND EIGENVECTORS OF A COMPLEX MATRIX
C     WHICH HAS BEEN REDUCED BY SUBROUTINE F01AMF TO UPPER
C     HESSENBERG FORM, H, STORED IN THE ARRAYS HR(N,N) AND
C     HI(N,N). THE REAL AND IMAGINARY PARTS OF THE EIGENVALUES
C     ARE FORMED IN THE ARRAYS WR(N) AND WI(N) RESPECTIVELY
C     AND THE UN-NORMALISED EIGENVECTORS ARE FORMED AS COLUMNS OF
C     THE ARRAYS VR(N,N) AND VI(N,N). LOW AND IUPP ARE TWO
C     INTEGERS PRODUCED IN BALANCING WHERE EIGENVALUES ARE
C     ISOLATED IN POSITIONS 1 TO LOW -1 AND IUPP+1 TO N. IF
C     BALANCING IS NOT USED LOW=1, IUPP=N. ACHEPS IS THE RELATIVE
C     MACHINE PRECISION. THE SUBROUTINE FAILS IF ALL EIGENVALUES
C     TAKE MORE THAN 30*N ITERATIONS.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02ARF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACHEPS
      INTEGER           IFAIL, IHI, IHR, IUPP, IVI, IVR, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  HI(IHI,N), HR(IHR,N), VI(IVI,N), VR(IVR,N),
     *                  WI(N), WR(N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AAHR, AHR, ANORM, SI, SR, TI, TR, XI, XR, XXI,
     *                  XXR, YI, YR, YYI, YYR, ZI, ZR
      INTEGER           I, I1, IEN, IEN1, II, IIEN, III, ISAVE, ITN,
     *                  ITS, IUPP1, IUPP11, J, J1, JJ, K, K1, KK, LOW1,
     *                  LOW11, M, M1, MM
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          A02AAF, A02ACF, F01AMZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
      ISAVE = IFAIL
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
      DO 40 I = 1, N
         DO 20 J = 1, N
            VR(I,J) = 0.0D0
            VI(I,J) = 0.0D0
   20    CONTINUE
         VR(I,I) = 1.0D0
   40 CONTINUE
      IUPP1 = IUPP - 1
      LOW1 = LOW + 1
      IF (LOW1.GT.IUPP1) GO TO 160
      DO 140 II = LOW1, IUPP1
         I = LOW1 + IUPP1 - II
         J = INTGER(I)
         I1 = I + 1
         IF (I1.GT.IUPP) GO TO 80
         DO 60 K = I1, IUPP
            III = I - 1
            VR(K,I) = HR(K,III)
            VI(K,I) = HI(K,III)
   60    CONTINUE
   80    IF (I.EQ.J) GO TO 140
         IF (I.GT.IUPP) GO TO 120
         DO 100 K = I, IUPP
            VR(I,K) = VR(J,K)
            VI(I,K) = VI(J,K)
            VR(J,K) = 0.0D0
            VI(J,K) = 0.0D0
  100    CONTINUE
  120    VR(J,I) = 1.0D0
  140 CONTINUE
  160 LOW11 = LOW - 1
      IF (1.GT.LOW11) GO TO 200
      DO 180 I = 1, LOW11
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  180 CONTINUE
  200 IUPP11 = IUPP + 1
      IF (IUPP11.GT.N) GO TO 240
      DO 220 I = IUPP11, N
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  220 CONTINUE
C     END OF ISOLATED ROOTS
  240 IEN = IUPP
  260 IF (IEN.LT.LOW) GO TO 800
      ITS = 0
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  280 IF (LOW1.GT.IEN) GO TO 320
      DO 300 KK = LOW1, IEN
         K = LOW1 + IEN - KK
         K1 = K - 1
         AHR = ABS(HR(K,K1)) + ABS(HI(K,K1))
         AAHR = ACHEPS*(ABS(HR(K1,K1))+ABS(HI(K1,K1))+ABS(HR(K,K))
     *          +ABS(HI(K,K)))
         IF (AHR.LE.AAHR) GO TO 340
  300 CONTINUE
  320 K = LOW
  340 IF (K.EQ.IEN) GO TO 780
      IF (ITN.LE.0) GO TO 1060
C     FORM SHIFT
      IF (ITS.EQ.10 .OR. ITS.EQ.20) GO TO 380
      SR = HR(IEN,IEN)
      SI = HI(IEN,IEN)
      IEN1 = IEN - 1
      XR = HR(IEN1,IEN)*HR(IEN,IEN1) - HI(IEN1,IEN)*HI(IEN,IEN1)
      XI = HR(IEN1,IEN)*HI(IEN,IEN1) + HI(IEN1,IEN)*HR(IEN,IEN1)
      IF (XR.EQ.0.0D0 .AND. XI.EQ.0.0D0) GO TO 400
      YR = (HR(IEN1,IEN1)-SR)/2.0D0
      YI = (HI(IEN1,IEN1)-SI)/2.0D0
      CALL A02AAF(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZR,ZI)
      IF ((YR*ZR+YI*ZI).GE.0.0D0) GO TO 360
      ZR = -ZR
      ZI = -ZI
  360 YYI = YI + ZI
      YYR = YR + ZR
      XXR = XR
      XXI = XI
      CALL A02ACF(XXR,XXI,YYR,YYI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 400
  380 IEN1 = IEN - 1
      SR = ABS(HR(IEN,IEN1)) + ABS(HR(IEN1,IEN-2))
      SI = ABS(HI(IEN,IEN1)) + ABS(HI(IEN1,IEN-2))
  400 IF (LOW.GT.IEN) GO TO 440
      DO 420 I = LOW, IEN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  420 CONTINUE
  440 TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
      J = K + 1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
      XR = ABS(HR(IEN1,IEN1)) + ABS(HI(IEN1,IEN1))
      YR = ABS(HR(IEN,IEN1)) + ABS(HI(IEN,IEN1))
      ZR = ABS(HR(IEN,IEN)) + ABS(HI(IEN,IEN))
      IF (J.GT.IEN1) GO TO 480
      DO 460 MM = J, IEN1
         M = J + IEN1 - MM
         YI = YR
         M1 = M - 1
         YR = ABS(HR(M,M1)) + ABS(HI(M,M1))
         XI = ZR
         ZR = XR
         XR = ABS(HR(M1,M1)) + ABS(HI(M1,M1))
         IF (YR.LE.(ACHEPS*ZR/YI*(ZR+XR+XI))) GO TO 500
  460 CONTINUE
  480 M = K
C     TRIANGULAR DECOMPOSITION H = L*R
  500 M1 = M + 1
      IF (M1.GT.IEN) GO TO 640
      DO 620 I = M1, IEN
         I1 = I - 1
         XR = HR(I1,I1)
         XI = HI(I1,I1)
         YR = HR(I,I1)
         YI = HI(I,I1)
         IF ((ABS(XR)+ABS(XI)).LT.(ABS(YR)+ABS(YI))) GO TO 520
         CALL A02ACF(YR,YI,XR,XI,ZR,ZI)
         WR(I) = -1.0D0
         GO TO 580
  520    IF (I1.GT.N) GO TO 560
C        INTERCHANGE ROWS OF HR AND HI
         DO 540 J = I1, N
            ZR = HR(I1,J)
            HR(I1,J) = HR(I,J)
            HR(I,J) = ZR
            ZI = HI(I1,J)
            HI(I1,J) = HI(I,J)
            HI(I,J) = ZI
  540    CONTINUE
  560    CALL A02ACF(XR,XI,YR,YI,ZR,ZI)
         WR(I) = 1.0D0
  580    HR(I,I1) = ZR
         HI(I,I1) = ZI
         IF (I.GT.N) GO TO 620
         DO 600 J = I, N
            HR(I,J) = HR(I,J) - ZR*HR(I1,J) + ZI*HI(I1,J)
            HI(I,J) = HI(I,J) - ZR*HI(I1,J) - ZI*HR(I1,J)
  600    CONTINUE
  620 CONTINUE
  640 IF (M1.GT.IEN) GO TO 280
C     COMPOSITION R*L = H
      DO 760 J = M1, IEN
         J1 = J - 1
         XR = HR(J,J1)
         XI = HI(J,J1)
         HR(J,J1) = 0.0D0
         HI(J,J1) = 0.0D0
C        INTERCHANGE COLUMNS OF HR, HI, VR, AND VI, IF NECESSARY
         IF (WR(J).LE.0.0D0) GO TO 700
         DO 660 I = 1, J
            ZR = HR(I,J1)
            HR(I,J1) = HR(I,J)
            HR(I,J) = ZR
            ZI = HI(I,J1)
            HI(I,J1) = HI(I,J)
            HI(I,J) = ZI
  660    CONTINUE
         IF (LOW.GT.IUPP) GO TO 700
         DO 680 I = LOW, IUPP
            ZR = VR(I,J1)
            VR(I,J1) = VR(I,J)
            VR(I,J) = ZR
            ZI = VI(I,J1)
            VI(I,J1) = VI(I,J)
            VI(I,J) = ZI
  680    CONTINUE
C        END INTERCHANGE COLUMNS
  700    DO 720 I = 1, J
            HR(I,J1) = HR(I,J1) + XR*HR(I,J) - XI*HI(I,J)
            HI(I,J1) = HI(I,J1) + XR*HI(I,J) + XI*HR(I,J)
  720    CONTINUE
         IF (LOW.GT.IUPP) GO TO 760
         DO 740 I = LOW, IUPP
            VR(I,J1) = VR(I,J1) + XR*VR(I,J) - XI*VI(I,J)
            VI(I,J1) = VI(I,J1) + XR*VI(I,J) + XI*VR(I,J)
  740    CONTINUE
C        END ACCUMULATE TRANSFORMATIONS
  760 CONTINUE
      GO TO 280
C     A ROOT FOUND
  780 WR(IEN) = HR(IEN,IEN) + TR
      WI(IEN) = HI(IEN,IEN) + TI
      IEN = IEN - 1
      GO TO 260
C     ALL ROOTS FOUND
  800 ANORM = 0.0D0
      DO 840 I = 1, N
         ANORM = ANORM + ABS(WR(I)) + ABS(WI(I))
         I1 = I + 1
         IF (I1.GT.N) GO TO 840
         DO 820 J = I1, N
            ANORM = ANORM + ABS(HR(I,J)) + ABS(HI(I,J))
  820    CONTINUE
  840 CONTINUE
      IF (ANORM.EQ.0.0D0) GO TO 1040
C     BACKSUBSTITUTE TO FIND VECTORS OF UPPER TRIANGULAR FORM
      IF (N.LT.2) GO TO 920
      DO 900 IIEN = 2, N
         IEN = 2 + N - IIEN
         XR = WR(IEN)
         XI = WI(IEN)
         IEN1 = IEN - 1
         DO 880 JJ = 1, IEN1
            J = 1 + IEN1 - JJ
            YR = XR - WR(J)
            YI = XI - WI(J)
            IF (YR.EQ.0.0D0 .AND. YI.EQ.0.0D0) YR = ACHEPS*ANORM
            CALL A02ACF(HR(J,IEN),HI(J,IEN),YR,YI,ZR,ZI)
            HR(J,IEN) = ZR
            HI(J,IEN) = ZI
            J1 = J - 1
            IF (J1.EQ.0) GO TO 880
            DO 860 I = 1, J1
               HR(I,IEN) = HR(I,IEN) + HR(I,J)*ZR - HI(I,J)*ZI
               HI(I,IEN) = HI(I,IEN) + HR(I,J)*ZI + HI(I,J)*ZR
  860       CONTINUE
  880    CONTINUE
  900 CONTINUE
C     END BACKSUB
C     MULTIPLY BY TRANSFORMATION MATRIX TO GIVE VECTORS OF
C     ORIGINAL FULL MATRIX
  920 IUPP1 = IUPP + 1
      DO 1000 J = 1, N
         M = MIN(J,LOW) - 1
         IF (M.LT.1) GO TO 960
         DO 940 I = 1, M
            VR(I,J) = HR(I,J)
            VI(I,J) = HI(I,J)
  940    CONTINUE
  960    M = J - 1
         IF (IUPP1.GT.M) GO TO 1000
         DO 980 I = IUPP1, M
            VR(I,J) = HR(I,J)
            VI(I,J) = HI(I,J)
  980    CONTINUE
 1000 CONTINUE
C     END VECTORS OF ISOLATED ROOTS
      IF (LOW.GT.N) GO TO 1040
      DO 1020 JJ = LOW, N
         J = LOW + N - JJ
         M = MIN(J-1,IUPP)
         IF (LOW.LE.M) CALL F01AMZ(VR(LOW,LOW),IVR,VI(LOW,LOW)
     *                             ,IVI,IUPP-LOW+1,M-LOW+1,HR(LOW,J)
     *                             ,1,HI(LOW,J),1,VR(LOW,J),VI(LOW,J))
 1020 CONTINUE
 1040 IFAIL = 0
      RETURN
 1060 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
