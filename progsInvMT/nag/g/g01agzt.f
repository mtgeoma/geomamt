      SUBROUTINE G01AGZ(ZMN,ZMX,NSTEP,ZNMIN,STEP,MAXA,MAXB)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     G01AGZ COMPUTES CHARACTERISTICS OF AN AXIS
C     FOR A GRAPH TO BE PLOTTED BY G01AGF
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STEP, ZMN, ZMX, ZNMIN
      INTEGER           MAXA, MAXB, NSTEP
C     .. Local Scalars ..
      DOUBLE PRECISION  AR, RINT, RNSTPZ, TENN, XMP, ZNM, ZNMAX
      INTEGER           I, II, IT, J, MINT
C     .. Local Arrays ..
      DOUBLE PRECISION  R(9)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, AINT, INT, LOG10, MAX, MIN, DBLE
C     .. Data statements ..
      DATA              R(1), R(2), R(3), R(4), R(5), R(6), R(7), R(8),
     *                  R(9)/0.1D0, 0.15D0, 0.2D0, 0.25D0, 0.4D0, 0.5D0,
     *                  0.6D0, 0.75D0, 0.8D0/
C     .. Executable Statements ..
      XMP = X02AJF()
      RNSTPZ = NSTEP
      RINT = (ZMX-ZMN)/(RNSTPZ+0.1D0)
   20 MINT = LOG10(RINT) - 2.0D0
      TENN = 10.0D0**MINT
      DO 60 J = 1, 11
         DO 40 I = 1, 9
            AR = R(I)
            IF (AR*TENN.GE.RINT) GO TO 80
   40    CONTINUE
         TENN = TENN*10.0D0
   60 CONTINUE
   80 STEP = TENN*AR
      IT = INT((1.0D0+2.0D0*XMP)*ZMN/STEP)
      AR = STEP*DBLE(IT)
  100 IF ((AR-STEP*0.05D0).LE.ZMN) GO TO 120
      AR = AR - STEP
      GO TO 100
  120 ZNMIN = AR
      ZNMAX = ZNMIN + STEP*(RNSTPZ+0.04D0)
      IF (ZNMAX.GE.ZMX) GO TO 140
      RINT = RINT*1.05D0
      GO TO 20
  140 MAXA = LOG10(ABS(ZNMAX)) + 1.0D0
      IF (ABS(ZNMIN).GT.ABS(ZNMAX)) MAXA = LOG10(ABS(ZNMIN)) + 1.0D0
      MAXA = MAX(0,MAXA)
      MAXB = LOG10(STEP) - 2.5D0
      MAXB = -MIN(MAXB,0)
      IF (MAXA+MAXB.GE.10) GO TO 180
      ZNM = ZNMIN
      DO 160 II = 1, 10
         I = II - 1
         IF (ZNM+STEP*(RNSTPZ+0.04D0).LT.ZMX) GO TO 180
         ZNMIN = ZNM
         ZNM = (ZNMIN*10.0D0**(MAXB-I))
         IF (ZNM.LT.0.0D0) ZNM = ZNM - 1.0D0
         ZNM = AINT(ZNM)/10.0D0**(MAXB-I)
  160 CONTINUE
  180 RETURN
      END
