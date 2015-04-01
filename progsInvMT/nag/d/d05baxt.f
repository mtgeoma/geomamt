      SUBROUTINE D05BAX(CK,CG,CF,ALIM,NXTPOL,TOL,THRESH,H,WT,STWT,LSTWT,
     *                  WKY1,WKY2,VF,VK,VG,WVG,NOUT,N,INFO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     ----------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      This subroutine computes the starting values to the required
C      accuracy. It uses either a 2-step ADAM or BDF to
C      evalute the solution with stepsizes H and H/2. This
C      procedure is continued when the required accuracy is
C      satsfied.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ----------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, H, THRESH, TOL
      INTEGER           INFO, LSTWT, N, NOUT, NXTPOL
C     .. Array Arguments ..
      DOUBLE PRECISION  STWT(LSTWT,0:1), VF(0:N), VG(0:2), VK(0:N),
     *                  WKY1(0:1), WKY2(0:2), WT(0:1), WVG(0:1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF, CG, CK
      EXTERNAL          CF, CG, CK
C     .. Scalars in Common ..
      DOUBLE PRECISION  GTOL, VGC
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST, DENMAX, GEVAL, GFUN, H2, H22, HALIM, HI,
     *                  HI2, SCALE, W0K0, Y0
      INTEGER           I, I2, IFLAG, IND, IR, NXNS
C     .. Local Arrays ..
      DOUBLE PRECISION  C(26)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          C05AXF, D05BAT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON            /AD05BA/GTOL
      COMMON            /BD05BA/VGC
C     .. Executable Statements ..
C
      HALIM = ALIM + H
      VK(0) = CK(0.D0)
      VK(1) = CK(H)
      VF(0) = CF(ALIM)
      VF(1) = CF(HALIM)
      W0K0 = WT(0)*VK(0)
C
C     ---- Evaluation of starting values with H  ----
C
      WKY1(0) = VF(0)
      VGC = CG(ALIM,WKY1(0))
      CALL D05BAT(CG,ALIM,H,STWT,LSTWT,WKY1,VF,VK,VG,INFO)
      IF (INFO.EQ.4) RETURN
C
C     ----  Start of extrapolation.  ----
C
      WKY2(0) = WKY1(0)
   20 CONTINUE
      WVG(0) = VG(0)
      WVG(1) = VG(1)
      NXTPOL = 2*NXTPOL
C
      IF (N.LT.NXTPOL*NOUT) THEN
         INFO = 5
         RETURN
      END IF
C
      NXNS = NXTPOL
      H2 = H/NXTPOL
      DO 40 I = NXNS/2, 1, -1
         I2 = 2*I
         VF(I2) = VF(I)
         VK(I2) = VK(I)
   40 CONTINUE
      HI2 = 2*H2
      HI = -H2
      DO 60 I = 1, NXNS, 2
         HI = HI + HI2
         VF(I) = CF(ALIM+HI)
         VK(I) = CK(HI)
   60 CONTINUE
C
C     ----  Evaluate the starting values with H/2**i.  ----
C
      CALL D05BAT(CG,ALIM,H2,STWT,LSTWT,WKY2,VF,VK,VG,INFO)
      IF (INFO.EQ.4) RETURN
C
      H22 = 2.D0*H2
      HALIM = H22 + ALIM
      SCALE = SQRT(X02AJF())
      IR = 1
      CONST = VF(2) + STWT(2,0)*VK(2)*VG(0) + STWT(2,1)*VK(1)*VG(1)
      IFLAG = 1
      Y0 = WKY2(1)
      IND = 1
   80 CALL C05AXF(Y0,GFUN,GTOL,IR,SCALE,C,IND,IFLAG)
      GEVAL = CG(HALIM,Y0)
      IF (IND.EQ.0) GO TO 100
      GFUN = Y0 - H2*W0K0*GEVAL - CONST
      GO TO 80
  100 CONTINUE
      IF (IFLAG.GT.0) THEN
         INFO = 4
         RETURN
      END IF
      WKY2(2) = Y0
      VG(2) = H2*GEVAL
      DENMAX = MAX(ABS(WKY1(1)),ABS(WKY2(2)),ABS(THRESH))
      WKY1(1) = (WKY2(2)-WKY1(1))/DENMAX
      IF (ABS(WKY1(1)).GT.TOL) THEN
         WKY1(1) = WKY2(1)
         GO TO 20
      END IF
      RETURN
      END
