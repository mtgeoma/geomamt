      SUBROUTINE D05BAV(CG,CF,CK,NXTPOL,VK,VG,WVG,VF,WKY1,WKY2,STWT,WT,
     *                  V,A,B,THETA,IQ,WWVG,WVK,LSTWT,L1WT,NSTART,H,IUP,
     *                  ALIM,TOL,THRESH,NOUT,N,INFO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine CARries on the EXTrapolation until the required
C     tolerance is satisfied.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, H, THRESH, TOL
      INTEGER           INFO, IQ, IUP, L1WT, LSTWT, N, NOUT, NSTART,
     *                  NXTPOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(6,0:5), B(0:IQ-1), STWT(LSTWT,0:NSTART),
     *                  THETA(0:IQ-1), V(0:NSTART), VF(0:N), VG(0:N),
     *                  VK(0:N), WKY1(0:N/2), WKY2(0:N), WT(0:L1WT),
     *                  WVG(0:N/2), WVK(0:NSTART-1,0:IQ-1,0:IQ-1),
     *                  WWVG(0:NSTART,0:IQ-1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF, CG, CK
      EXTERNAL          CF, CG, CK
C     .. Scalars in Common ..
      DOUBLE PRECISION  VGC
C     .. Local Scalars ..
      DOUBLE PRECISION  H2, HI, HI2, VG0
      INTEGER           I, I2, INDEX, IUPD2
C     .. External Subroutines ..
      EXTERNAL          D05BAK, D05BAP, D05BAQ, D05BAR, D05BAS, D05BAT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /BD05BA/VGC
C     .. Executable Statements ..
C
   20 CONTINUE
      HI = H/NXTPOL
      NXTPOL = 2*NXTPOL
C
      IF (N.LT.NXTPOL*NOUT) THEN
         INFO = 5
         RETURN
      END IF
C
      IUPD2 = IUP
      DO 40 I = 0, IUPD2
         WKY1(I) = WKY2(I)
         WVG(I) = VG(I)
   40 CONTINUE
      IUP = 2*IUP
      H2 = H/NXTPOL
      VG0 = H2*VGC
      VG(0) = VG0
      DO 60 I = IUPD2, 1, -1
         I2 = 2*I
         VF(I2) = VF(I)
         VK(I2) = VK(I)
   60 CONTINUE
      HI2 = 2*H2
      HI = -H2
      DO 80 I = 1, IUP, 2
         HI = HI + HI2
         VF(I) = CF(ALIM+HI)
         VK(I) = CK(HI)
   80 CONTINUE
C
C     ---- Evaluate the starting values with H/2**i. ----
C
      IF (NSTART.EQ.1) THEN
         CALL D05BAT(CG,ALIM,H2,STWT,LSTWT,WKY2,VF,VK,VG,INFO)
         IF (INFO.EQ.4) RETURN
      ELSE
         CALL D05BAK(CK,H2,THETA,VK,WVK,IQ,NSTART)
C
         CALL D05BAP(CG,CF,ALIM,WVK,WWVG,WKY2,VF,VG,H2,THETA,B,A,IQ,
     *               NSTART)
      END IF
C
C     ----  Evaluate the inhomogenious trems. ----
C
      CALL D05BAQ(WKY2,STWT,VK,VG,VF,V,LSTWT,NSTART,NSTART+1,IUP,1)
C
C     ----  Solve the convolution equation. ----
C
      CALL D05BAS(CG,H2,ALIM,WT,L1WT,VK,VG,WKY2,NSTART,IUP,NSTART,1,
     *            INFO)
      IF (INFO.EQ.4) RETURN
C
C     ---- Check the errors  ----
C
      CALL D05BAR(WKY1(1),WKY2(1),THRESH,IUPD2,INDEX)
      IF (ABS(WKY1(INDEX)).GT.TOL) GO TO 20
C
      RETURN
      END
