      SUBROUTINE D02GBF(A,B,N,TOL,FCNF,FCNG,C,D,GAM,MNP,X,Y,NP,W,LW,IW,
     *                  LIW,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-305 (SEP 1981).
C     MARK 9C REVISED. IER-373 (JUN 1982)
C     MARK 10 REVISED. IER-377 (JUN 1982).
C     MARK 10A REVISED. IER-389 (OCT 1982).
C     MARK 10B REVISED. IER-400 (JAN 1983).
C     MARK 11 REVISED. IER-434 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. IER-606 (APR 1988).
C     FCNF, FCNG
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02GBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, TOL
      INTEGER           IFAIL, LIW, LW, MNP, N, NP
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N,N), D(N,N), GAM(N), W(LW), X(MNP), Y(N,MNP)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          FCNF, FCNG
C     .. Local Scalars ..
      DOUBLE PRECISION  DELEPS, TEMP
      INTEGER           I, I1, IFAIL1, IND, IND1, INIT, J, JJ, K, LIN,
     *                  LP, LWORK, M1, M2, M3, MP, N1, NI, NSQ,
     *                  NUMBEG, NUMMIX
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02GAX, D02GAY, D02GAZ, D02GBZ, D02HBS, D02RAZ,
     *                  M01DBF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      NSQ = N*N
      IF (B.GT.A .AND. (NP.EQ.0 .OR. NP.GE.4)
     *    .AND. TOL.GT.0.0D0 .AND. MNP.GE.32 .AND. N.GT.0 .AND. NP.LE.
     *    MNP .AND. LW.GE.MNP*(3*NSQ+5*N+2)
     *    +3*NSQ+5*N .AND. LIW.GE.MNP*(2*N+1)+N) GO TO 40
C     INPUT ERROR
   20 IND = 1
      GO TO 420
   40 DO 80 I = 1, N
         K = 0
         DO 60 J = 1, N
            IF (C(I,J).NE.0.0D0 .AND. K.EQ.3) K = 2
            IF (C(I,J).NE.0.0D0 .AND. K.EQ.0) K = 1
            IF (D(I,J).NE.0.0D0 .AND. K.EQ.1) K = 2
            IF (D(I,J).NE.0.0D0 .AND. K.EQ.0) K = 3
   60    CONTINUE
         IW(I) = K
         IF (IW(I).NE.0) GO TO 80
         IND = 2
         IW(1) = I
         GO TO 420
   80 CONTINUE
      JJ = 0
      DO 120 I = 1, N
         NI = N + I
         W(I) = 0.D0
         W(NI) = 0.D0
         DO 100 J = 1, N
            IF (C(J,I).NE.0.D0) W(I) = 1.D0
            IF (D(J,I).NE.0.D0) W(NI) = 1.D0
  100    CONTINUE
         IF (W(I).NE.0.D0) JJ = JJ + 1
         IF (W(NI).NE.0.D0) JJ = JJ + 1
  120 CONTINUE
      IF (JJ.GE.N) GO TO 140
      IND = 2
      IW(1) = -JJ
      GO TO 420
  140 CALL D02GBZ(C,D,N,W,W(NSQ+1),IND1,W(2*NSQ+1))
      IF (IND1.EQ.0) GO TO 160
      IW(1) = IND1
      IND = 5
      GO TO 420
  160 NUMBEG = 0
      NUMMIX = 0
      DO 180 I = 1, N
         IF (IW(I).EQ.1) NUMBEG = NUMBEG + 1
         IF (IW(I).EQ.2) NUMMIX = NUMMIX + 1
  180 CONTINUE
      IF (NUMMIX.EQ.N) GO TO 300
      I1 = NUMBEG + NUMMIX
      IF (NUMBEG.NE.N .AND. I1.NE.0) GO TO 200
      IND = 2
      IW(1) = 0
      GO TO 420
  200 IFAIL1 = 1
      N1 = N + 1
      CALL M01DBF(IW,1,N,'A',IW(N1),IFAIL1)
      DO 280 I = 1, N
         NI = N + I
         DO 220 K = I, N
            J = N + K
            IF (IW(J).EQ.I) GO TO 240
  220    CONTINUE
  240    CONTINUE
         IF (K.EQ.I) GO TO 280
         DO 260 JJ = 1, N
            TEMP = C(K,JJ)
            C(K,JJ) = C(I,JJ)
            C(I,JJ) = TEMP
            TEMP = D(K,JJ)
            D(K,JJ) = D(I,JJ)
            D(I,JJ) = TEMP
  260    CONTINUE
         TEMP = GAM(K)
         GAM(K) = GAM(I)
         GAM(I) = TEMP
         IW(J) = IW(NI)
  280 CONTINUE
  300 IF (NP.EQ.0) GO TO 340
      IF (A.NE.X(1) .OR. B.NE.X(NP)) GO TO 20
      N1 = NP - 1
      DO 320 I = 1, N1
         IF (X(I+1).LE.X(I)) GO TO 20
  320 CONTINUE
      INIT = 1
      GO TO 360
  340 INIT = 0
      NP = 4
  360 LIN = 1
      DELEPS = 0.0D0
      LP = MOD(IFAIL/10,10)
      MP = MOD(IFAIL/100,10)
      M1 = 1 + N
      M2 = M1 + N*N
      M3 = M2 + N
      LWORK = LW - M3 + 1
      CALL D02RAZ(N,MNP,NP,NUMBEG,NUMMIX,A,B,TOL,X,Y,N,W(1)
     *            ,D02HBS,D02GAZ,D02GAX,FCNF,FCNG,D02GAZ,D02GAY,D02GAZ,
     *            D02GAX,W(M1),W(M2),C,D,GAM,W,W,W(M3)
     *            ,LWORK,IW,LIW,DELEPS,LP,MP,INIT,LIN,IND)
      IF (IND.EQ.0) GO TO 420
      GO TO (380,400,400,400,380,380,380,380,380) IND
  380 IND = 4
      GO TO 420
  400 IND = 3
  420 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
