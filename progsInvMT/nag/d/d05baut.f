      SUBROUTINE D05BAU(CK,CG,CF,ALIM,NSTART,L1WT,LSTWT,NXTPOL,TOL,
     *                  THRESH,H,WT,STWT,V,WKY1,WKY2,VF,VK,VG,WVG,A,B,
     *                  THETA,IQ,WWVG,WVK,NOUT,N,INFO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     ----------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      This subroutine computes the starting values to the required
C      accuracy. With a stepsize H, it uses a Runge-Kutte scheme to
C      evaluate the starting values an uses a reducible LMM to
C      evalute the remaining solution with stepsize H/2. This
C      procedure is continued when the required accuracy is
C      satsfied.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, H, THRESH, TOL
      INTEGER           INFO, IQ, L1WT, LSTWT, N, NOUT, NSTART, NXTPOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(6,0:5), B(0:IQ-1), STWT(LSTWT,0:NSTART),
     *                  THETA(0:IQ-1), V(0:NSTART), VF(0:N),
     *                  VG(0:2*NSTART), VK(0:N), WKY1(0:NSTART),
     *                  WKY2(0:2*NSTART), WT(0:L1WT), WVG(0:NSTART),
     *                  WVK(0:NSTART-1,0:IQ-1,0:IQ-1),
     *                  WWVG(0:NSTART,0:IQ-1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF, CG, CK
      EXTERNAL          CF, CG, CK
C     .. Scalars in Common ..
      DOUBLE PRECISION  VGC
C     .. Local Scalars ..
      DOUBLE PRECISION  H2, HI, HI2, VG0
      INTEGER           I, I2, INDEX, NSTRT2, NXNS
C     .. External Subroutines ..
      EXTERNAL          D05BAK, D05BAP, D05BAQ, D05BAR, D05BAS
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /BD05BA/VGC
C     .. Executable Statements ..
C
      NSTRT2 = 2*NSTART
      HI = 0.D0
      VK(0) = CK(HI)
      DO 20 I = 1, NSTART
         HI = HI + H
         VK(I) = CK(HI)
   20 CONTINUE
C
C     ---- Evaluation of starting values with H  ----
C
      HI = ALIM
      VF(0) = CF(HI)
      DO 40 I = 1, NSTART
         HI = HI + H
         VF(I) = CF(HI)
   40 CONTINUE
      WKY1(0) = VF(0)
      VGC = CG(ALIM,WKY1(0))
      VG0 = H*VGC
      VG(0) = VG0
      CALL D05BAK(CK,H,THETA,VK,WVK,IQ,NSTART)
C
      CALL D05BAP(CG,CF,ALIM,WVK,WWVG,WKY1,VF,VG,H,THETA,B,A,IQ,NSTART)
C
C     ----  Start of extrapolation.  ----
C
      WKY2(0) = WKY1(0)
   60 CONTINUE
      DO 80 I = 0, NSTART
         WVG(I) = VG(I)
   80 CONTINUE
      NXTPOL = 2*NXTPOL
C
      IF (N.LT.NXTPOL*NOUT) THEN
         INFO = 5
         RETURN
      END IF
C
      NXNS = NSTART*NXTPOL
      H2 = H/NXTPOL
      VG0 = H2*VGC
      VG(0) = VG0
      DO 100 I = NXNS/2, 1, -1
         I2 = 2*I
         VF(I2) = VF(I)
         VK(I2) = VK(I)
  100 CONTINUE
      HI2 = 2*H2
      HI = -H2
      DO 120 I = 1, NXNS, 2
         HI = HI + HI2
         VF(I) = CF(ALIM+HI)
         VK(I) = CK(HI)
  120 CONTINUE
C
C     ----  Evaluate the starting values with H/2**i.  ----
C
      CALL D05BAK(CK,H2,THETA,VK,WVK,IQ,NSTART)
C
      CALL D05BAP(CG,CF,ALIM,WVK,WWVG,WKY2,VF,VG,H2,THETA,B,A,IQ,NSTART)
C
C     ----  Evaluate the inhomogenious trems.  ----
C
      CALL D05BAQ(WKY2,STWT,VK,VG,VF,V,LSTWT,NSTART,NSTART+1,NSTRT2,1)
C
C     ----  Solve the convolution equation. ----
C
      CALL D05BAS(CG,H2,ALIM,WT,L1WT,VK,VG,WKY2,NSTART,NSTRT2,NSTART,1,
     *            INFO)
      IF (INFO.EQ.4) RETURN
C
C     ---- Check the errors ----
C
      CALL D05BAR(WKY1(1),WKY2(1),THRESH,NSTART,INDEX)
      IF (ABS(WKY1(INDEX)).GT.TOL) THEN
         DO 140 I = 1, NSTART
            WKY1(I) = WKY2(I)
  140    CONTINUE
         GO TO 60
      END IF
      RETURN
      END
