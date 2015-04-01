      SUBROUTINE G08BAF(X,N,N1,R,ITEST,W,V,PW,PV,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 10B REVISED. IER-404 (JAN 1983).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G08DAF PERFORMS MOODS AND DAVIDS TESTS FOR DISPERSION
C     DIFFERENCES BETWEEN TWO INDEPENDENT SAMPLES OF
C     POSSIBLY UNEQUAL SIZE.
C
C     REFERENCE  - B.E.COOPER - STATISTICS FOR EXPERIMENTALISTS
C     AUTHOR     - J. LLOYD-JONES (U.M.R.C.C.)
C
C     PARAMETERS -
C     X     - DATA FOR TWO SAMPLES
C     N     - TOTAL SIZE OF BOTH SAMPLES
C     N1    - SIZE OF FIRST SAMPLE
C     R     - VECTOR OF RANKS
C     ITEST - TEST REQUIRED - SET TO 1 FOR DAVIDS TEST
C                                    2 FOR MOODS TEST
C                                    0 FOR BOTH
C     W     - MOODS MEASURE OF DISPERSION DIFFERENCES
C     V     - DAVIDS MEASURE OF DISPERSION DIFFERENCES
C     PW    - SIGNIFICANCE OF W
C     PV    - SIGNIFICANCE OF V
C
C     CHECK PARAMETERS
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PV, PW, V, W
      INTEGER           IFAIL, ITEST, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EV, EW, FN, FN1, FN2, VARV, VARW, VV, VX, WW,
     *                  WX, XF, XMD, XMM
      INTEGER           I, IERROR, IFAWK, N2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF
      INTEGER           P01ABF
      EXTERNAL          S15ABF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08AEZ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IF (N.LE.2) GO TO 120
      IF (N1.LE.1 .OR. N1.GE.N) GO TO 140
      IF (ITEST.LT.0 .OR. ITEST.GT.2) GO TO 160
C     NO PARAMETER ERRORS
      FN = N
      FN1 = N1
      N2 = N - N1
      FN2 = N2
C     USE G08AEZ TO FIND RANKS OF POOLED SAMPLE
C     XF IS A DUMMY PARAMETER OF G08AEZ AND IS USED NOWHERE ELSE
      CALL G08AEZ(X,R,N,0,XF)
C     ********************************************
C     CALCULATE TEST STATISTIC(S) AND SIGNIFICANCE
C     ********************************************
      IF (ITEST.EQ.1) GO TO 40
C     PERFORM MOODS TEST
      XMM = (FN+1.0D0)*0.5D0
      W = 0.0D0
      DO 20 I = 1, N1
         WW = R(I) - XMM
         W = W + WW*WW
   20 CONTINUE
      EW = FN1*(FN*FN-1.0D0)/12.0D0
      VARW = FN1*FN2*(FN+1.0D0)*(FN*FN-4.0D0)/180.0D0
      WX = (W-EW)/(SQRT(VARW))
      IFAWK = 1
      PW = S15ABF(WX,IFAWK)
   40 IF (ITEST.EQ.2) GO TO 100
C     PERFORM DAVIDS TEST
      XMD = 0.0D0
      DO 60 I = 1, N1
         XMD = XMD + R(I)
   60 CONTINUE
      XMD = XMD/FN1
      V = 0.0D0
      DO 80 I = 1, N1
         VV = R(I) - XMD
         V = V + VV*VV
   80 CONTINUE
      V = V/(FN1-1.0D0)
      EV = FN*(FN+1.0D0)/12.0D0
      VARV = FN2*FN*(FN+1.0D0)*(3.0D0*(FN+1.0D0)*(FN1+1.0D0)-FN*FN1)
     *       /(360.0D0*FN1*(FN1-1.0D0))
      VX = (V-EV)/(SQRT(VARV))
      IFAWK = 1
      PV = S15ABF(VX,IFAWK)
  100 IFAIL = 0
      GO TO 200
C     ********
C     ERRORS
C     ********
  120 IERROR = 1
      GO TO 180
  140 IERROR = 2
      GO TO 180
  160 IERROR = 3
  180 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
  200 RETURN
      END
