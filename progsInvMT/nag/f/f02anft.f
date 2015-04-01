      SUBROUTINE F02ANF(NN,ACC,HR,IHR,HI,IHI,WR,WI,LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 9 REVISED. IER-326 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMLR
C     FINDS THE EIGENVALUES OF A COMPLEX UPPER-HESSENBERG MATRIX,
C     H, STORED IN THE ARRAYS HR(N,N) AND HI(N,N), AND
C     STORES THE REAL PARTS IN THE ARRAY WR(N) AND THE COMPLEX
C     PARTS IN THE ARRAY WI(N). ACC IS THE RELATIVE MACHINE
C     PRECISION. THE SUBROUTINE FAILS IF ALL EIGENVALUES TAKE MORE
C     THAN 30*N ITERATIONS.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02ANF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC
      INTEGER           IHI, IHR, LFAIL, NN
C     .. Array Arguments ..
      DOUBLE PRECISION  HI(IHI,NN), HR(IHR,NN), WI(NN), WR(NN)
C     .. Local Scalars ..
      DOUBLE PRECISION  SI, SR, TI, TR, XI, XR, YI, YR, ZI, ZR, ZXI, ZXR
      INTEGER           I, I1, ISAVE, ITN, ITS, J, L, LL, M, M1, MM, N,
     *                  N1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          A02AAF, A02ACF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 0
      N = NN
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
   20 IF (N.EQ.0) GO TO 480
      ITS = 0
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
   40 L = N + 1
      IF (N.EQ.1) GO TO 80
      DO 60 LL = 2, N
         L = L - 1
         IF ((ABS(HR(L,L-1))+ABS(HI(L,L-1))).LE.(ACC*(ABS(HR(L-1,L-1))
     *       +ABS(HI(L-1,L-1))+ABS(HR(L,L))+ABS(HI(L,L))))) GO TO 100
   60 CONTINUE
   80 L = 1
  100 IF (L.EQ.N) GO TO 460
      IF (ITN.GT.0) GO TO 120
      LFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
C     FORM SHIFT
  120 IF ((ITS.NE.10) .AND. (ITS.NE.20)) GO TO 140
      SR = ABS(HR(N,N-1)) + ABS(HR(N-1,N-2))
      SI = ABS(HI(N,N-1)) + ABS(HI(N-1,N-2))
      GO TO 180
  140 SR = HR(N,N)
      SI = HI(N,N)
      N1 = N - 1
      XR = HR(N1,N)*HR(N,N1) - HI(N1,N)*HI(N,N1)
      XI = HR(N1,N)*HI(N,N1) + HI(N1,N)*HR(N,N1)
      IF ((XR.EQ.0.0D0) .AND. (XI.EQ.0.0D0)) GO TO 180
      YR = (HR(N-1,N-1)-SR)/2.0D0
      YI = (HI(N-1,N-1)-SI)/2.0D0
      CALL A02AAF(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZR,ZI)
      IF ((YR*ZR+YI*ZI).GE.0.0D0) GO TO 160
      ZR = -ZR
      ZI = -ZI
  160 CALL A02ACF(XR,XI,YR+ZR,YI+ZI,ZXR,ZXI)
      XR = ZXR
      XI = ZXI
      SR = SR - XR
      SI = SI - XI
  180 DO 200 I = 1, N
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  200 CONTINUE
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
      J = L + 1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
      XR = ABS(HR(N-1,N-1)) + ABS(HI(N-1,N-1))
      YR = ABS(HR(N,N-1)) + ABS(HI(N,N-1))
      ZR = ABS(HR(N,N)) + ABS(HI(N,N))
      M = N
      N1 = N - 1
      IF (J.GT.N1) GO TO 240
      DO 220 MM = J, N1
         M = M - 1
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZR
         ZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         IF (YR.LE.(ACC*ZR/YI*(ZR+XR+XI))) GO TO 260
  220 CONTINUE
  240 M = L
C     TRIANGULAR DECOMPOSITION H = L*R
  260 M1 = M + 1
      DO 360 I = M1, N
         XR = HR(I-1,I-1)
         XI = HI(I-1,I-1)
         YR = HR(I,I-1)
         YI = HI(I,I-1)
         IF ((ABS(XR)+ABS(XI)).GE.(ABS(YR)+ABS(YI))) GO TO 300
C        INTERCHANGE ROWS OF HR AND HI
         I1 = I - 1
         DO 280 J = I1, N
            ZR = HR(I-1,J)
            HR(I-1,J) = HR(I,J)
            HR(I,J) = ZR
            ZI = HI(I-1,J)
            HI(I-1,J) = HI(I,J)
            HI(I,J) = ZI
  280    CONTINUE
         CALL A02ACF(XR,XI,YR,YI,ZR,ZI)
         WR(I) = 1.0D0
         GO TO 320
  300    CALL A02ACF(YR,YI,XR,XI,ZR,ZI)
         WR(I) = -1.0D0
  320    HR(I,I-1) = ZR
         HI(I,I-1) = ZI
         DO 340 J = I, N
            HR(I,J) = HR(I,J) - ZR*HR(I-1,J) + ZI*HI(I-1,J)
            HI(I,J) = HI(I,J) - ZR*HI(I-1,J) - ZI*HR(I-1,J)
  340    CONTINUE
  360 CONTINUE
C     COMPOSITION R*L = H
      M1 = M + 1
      DO 440 J = M1, N
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0D0
         HI(J,J-1) = 0.0D0
C        INTERCHANGE COLUMNS OF HR AND HI, IF NECESSARY
         IF (WR(J).LE.0.0D0) GO TO 400
         DO 380 I = L, J
            ZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZR
            ZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZI
  380    CONTINUE
  400    DO 420 I = L, J
            HR(I,J-1) = HR(I,J-1) + XR*HR(I,J) - XI*HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR*HI(I,J) + XI*HR(I,J)
  420    CONTINUE
  440 CONTINUE
      GO TO 40
C     A ROOT FOUND
  460 WR(N) = HR(N,N) + TR
      WI(N) = HI(N,N) + TI
      N = N - 1
      GO TO 20
  480 RETURN
      END
