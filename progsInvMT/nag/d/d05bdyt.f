      SUBROUTINE D05BDY(KC,FC,GC,YS,N,FFT,LENP,N1,S,ISM1,H,TOLNL,WT,SWT,
     *                  VKC,VG,FTVK,WKS,WKVG,FJAC,A,WN,NCT,FORS,IFLAG)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine computes  the position of blocks of lags which
C     is required for computation of the lag using FFT techniques.
C     It also computes the solution. The DFT of the kernel is
C     evaluated once and stored for subsequent use.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, TOLNL
      INTEGER           IFLAG, ISM1, LENP, N, N1, S
      CHARACTER         FFT, FORS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(ISM1,ISM1), FJAC(ISM1,ISM1), FTVK(2*N1),
     *                  SWT(N1+S,0:ISM1), VG(0:N1+ISM1), VKC(N-1),
     *                  WKS(N1+2*ISM1), WKVG(N1), WN(ISM1,4), WT(N1),
     *                  YS(0:N1+ISM1)
      INTEGER           NCT(*)
C     .. Function Arguments ..
      DOUBLE PRECISION  FC, GC, KC
      EXTERNAL          FC, GC, KC
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSFUN, FACTOR, HI, SQRH, SUM, VG0, WK0, YTOL
      INTEGER           I, I1, IFAIL, IFAIL1, IIPOS, IJ, IL, ILN, INTS,
     *                  IPOS, IPS1, IQ, IREVCM, ISEC, IU, IUL, IUL2,
     *                  IUN, IW, J, JJ, LDFJAC, LR, M, ML, MODE, MU,
     *                  NJJ, NJJ2, NNBLK, NOBLK
C     .. Local Arrays ..
      DOUBLE PRECISION  DIAG(11), FUN(11), QTF(11), R(66), TS(11), V(11)
C     .. External Subroutines ..
      EXTERNAL          C05NDF, C06EAF, D05BDX, D05BYW
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, LOG, NINT, SQRT
C     .. Executable Statements ..
C
      SQRH = SQRT(H)
C
C     ... Evaluation of the kernel and the function  f. ...
C
      DO 20 I = 1, S - 2
         WKS(I) = KC(DBLE(I-ISM1)*H)
   20 CONTINUE
C
      IF (FORS.EQ.'F') GO TO 40
C
      YS(0) = FC(0.D0)
   40 WKS(ISM1) = KC(0.D0)
C
      DO 60 I = 1, N1 + ISM1
         HI = H*DBLE(I)
         WKS(I+ISM1) = KC(HI)
         YS(I) = FC(HI)
   60 CONTINUE
C
C     ... Evaluation of starting values. ...
C
      VG(0) = SQRH*GC(0.0D0,YS(0))
      VG0 = VG(0)
C
C     ... Solving the nonlinear system for starting values . ...
C
      DO 100 I = 1, ISM1
C
         I1 = I + 1
         IPS1 = I + ISM1
         V(I) = VG0*SWT(I1,0)*WKS(IPS1) + YS(I)
         TS(I) = DBLE(I)*H
C
         DO 80 J = 1, ISM1
            A(I,J) = SQRH*SWT(I1,J)*WKS(IPS1-J)
   80    CONTINUE
C
  100 CONTINUE
C
      DO 120 I = 1, ISM1
         DIAG(I) = 1.D0
  120 CONTINUE
      ML = ISM1 - 1
      MU = ISM1 - 1
      EPSFUN = 0.D0
      FACTOR = 100.D0
      MODE = 2
      YTOL = TOLNL
      LDFJAC = ISM1
      LR = ISM1*(ISM1+1)/2
      IFAIL1 = 1
      IREVCM = 0
C
  140 CALL C05NDF(IREVCM,ISM1,YS(1),FUN,YTOL,ML,MU,EPSFUN,DIAG,MODE,
     *            FACTOR,FJAC,LDFJAC,R,LR,QTF,WN,IFAIL1)
C
      IF (IREVCM.EQ.1) THEN
         GO TO 140
      ELSE IF (IREVCM.EQ.2) THEN
C
         DO 160 I = 1, ISM1
            VG(I) = GC(TS(I),YS(I))
  160    CONTINUE
C
         DO 200 I = 1, ISM1
C
            IF (FORS.EQ.'S') THEN
               SUM = YS(I) - V(I)
            ELSE IF (FORS.EQ.'F') THEN
               SUM = -V(I)
            END IF
C
            DO 180 J = 1, ISM1
               SUM = SUM - A(I,J)*VG(J)
  180       CONTINUE
C
            FUN(I) = SUM
C
  200    CONTINUE
C
         GO TO 140
C
      ELSE
         IF (IFAIL1.NE.0) THEN
            IFLAG = 2
            RETURN
         END IF
      END IF
C
C     ... Evaluation of inhomogenious terms. ...
C
      DO 220 I = 1, ISM1
         VG(I) = SQRH*GC(TS(I),YS(I))
  220 CONTINUE
C
      DO 260 I = S, N1 + ISM1
         SUM = 0.0D0
         I1 = I + 1
C
         DO 240 J = 0, ISM1
            SUM = SUM + SWT(I1,J)*WKS(I-J+ISM1)*VG(J)
  240    CONTINUE
C
         YS(I) = YS(I) + SUM
  260 CONTINUE
C
C     ... Absorbing the convolution weights ...
C     ... in the values of kernel.          ...
C
      DO 280 I = 1, N1
         WKS(I) = SQRH*WKS(I+ISM1-1)*WT(I)
  280 CONTINUE
C
      WK0 = WKS(1)
C
      DO 300 I = 1, N - 1
         VKC(I) = WKS(I+1)
  300 CONTINUE
C
      IF (FFT.EQ.'N') THEN
         INTS = 0
         CALL D05BDX(GC,YS,N,INTS,N1,S,H,VG,VKC,WK0,FORS,TOLNL,IFLAG)
         GO TO 520
      END IF
C
C     ... Evaluate the DFT of the kernel once. ...
C
      IQ = NINT(LOG(DBLE(N))/LOG(DBLE(2)))
      NOBLK = LENP + 1 - IQ
      IPOS = 0
      NNBLK = 0
      IFAIL = 0
C
      DO 340 M = 1, NOBLK
         IPOS = NNBLK + IPOS
         NNBLK = N*(2**M)
C
         DO 320 I = 1, NNBLK
            FTVK(IPOS+I) = WKS(I)
  320    CONTINUE
C
         CALL C06EAF(FTVK(IPOS+1),NNBLK,IFAIL)
  340 CONTINUE
C
C     ... Start the solution of  Yj,  j = s,...,N + s -1. ...
C
      INTS = 0
C
      CALL D05BDX(GC,YS,N,INTS,N1,S,H,VG,VKC,WK0,FORS,TOLNL,IFLAG)
C
      IF (IFLAG.NE.0) GO TO 520
C
      NJJ2 = 1
      IPOS = 0
C
C     ... Start of the loop for evaluation of the ...
C     ... lag in big blocks by FFT.                      ...
C
      DO 500 M = 1, NOBLK
         IPOS = NJJ2 + IPOS
         JJ = 2**(M-1)
         NJJ = N*JJ
         NJJ2 = NJJ*2
         NCT(JJ) = 0
C
         DO 360 I = 1, NJJ
            WKVG(I) = VG(I+ISM1)
  360    CONTINUE
C
         DO 380 I = NJJ + 1, NJJ2
            WKVG(I) = 0.0D0
  380    CONTINUE
C
         ISEC = 1
         CALL D05BYW(WKVG,FTVK(IPOS),NJJ2,ISEC)
C
         DO 400 I = NJJ + S, NJJ2 + ISM1
            YS(I) = YS(I) + WKVG(I-ISM1)
  400    CONTINUE
C
         CALL D05BDX(GC,YS,N,NJJ,N1,S,H,VG,VKC,WK0,FORS,TOLNL,IFLAG)
C
         IF (IFLAG.NE.0) GO TO 520
C
C        ... Start of loop for evaluation of lag in ...
C        ... subsequent small blocks.               ...
C
         DO 480 IJ = 1, JJ - 1
            IW = NCT(IJ)
            IL = JJ + IW
            IU = JJ + IJ
            ILN = IL*N
            IUN = IU*N
            IUL = IUN - ILN
            IUL2 = IUL*2
            NCT(IU) = IL
C
            DO 420 I = 1, IUL
               WKVG(I) = VG(I+ILN+ISM1)
  420       CONTINUE
C
            DO 440 I = IUL + 1, IUL2
               WKVG(I) = 0.0D0
  440       CONTINUE
C
            IF (IUL.EQ.N) THEN
               IIPOS = 1
            ELSE
               IIPOS = IUL2 - 2*N + 1
            END IF
C
            ISEC = 1
            CALL D05BYW(WKVG,FTVK(IIPOS),IUL2,ISEC)
C
            DO 460 I = IUN + S, IUN + IUL + ISM1
               YS(I) = YS(I) + WKVG(I-ILN-ISM1)
  460       CONTINUE
C
            CALL D05BDX(GC,YS,N,IUN,N1,S,H,VG,VKC,WK0,FORS,TOLNL,IFLAG)
C
            IF (IFLAG.NE.0) GO TO 520
C
  480    CONTINUE
C
  500 CONTINUE
  520 CONTINUE
      RETURN
      END
