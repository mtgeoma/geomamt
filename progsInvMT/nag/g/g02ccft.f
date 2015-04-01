      SUBROUTINE G02CCF(N,X,Y,XMISS,YMISS,RESULT,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 8 REVISED. IER-229 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     NAG SUBROUTINE G02CCF
C     WRITTEN  6.10.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     PERFORMS A SIMPLE LINEAR REGRESSION WITH DEPENDENT VARIABLE Y
C     AND INDEPENDENT VARIABLE X, OMITTING CASES INVOLVING MISSING
C     VALUES.
C
C     USES NAG ERROR ROUTINE P01AAF
C     NAG LIBRARY ROUTINES   X02AKF, X02ALF
C     X02BEF
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE-DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02CCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMISS, YMISS
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RESULT(21), X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, BIGNUM, FN, R, RELXM, RELYM, SMALL, STDX,
     *                  STDY, XBAR, XI, YBAR, YI
      INTEGER           I, IERROR
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AKF, X02ALF
      INTEGER           P01ABF, X02BEF
      EXTERNAL          X02AKF, X02ALF, P01ABF, X02BEF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      ACC = 0.1D0**(X02BEF(ACC)-2)
      BIGNUM = X02ALF()
      SMALL = X02AKF()
      IERROR = 1
      IF (N.LT.3) GO TO 240
      RELXM = ABS(ACC*XMISS)
      RELYM = ABS(ACC*YMISS)
      FN = 0.0D0
      XBAR = 0.0D0
      YBAR = 0.0D0
      DO 20 I = 1, N
         IF (ABS(X(I)-XMISS).LE.RELXM .OR. ABS(Y(I)-YMISS).LE.RELYM)
     *       GO TO 20
         FN = FN + 1.0D0
         XBAR = XBAR + X(I)
         YBAR = YBAR + Y(I)
   20 CONTINUE
C
C     WANT TO TEST FN.LT.3.0, BUT TO ALLOW FOR MACHINE ROUNDING
C     WE ACTUALLY TEST FN.LT.2.5
C
      IF (FN.LT.2.5D0) GO TO 200
      XBAR = XBAR/FN
      YBAR = YBAR/FN
      STDX = 0.0D0
      STDY = 0.0D0
      R = 0.0D0
      DO 40 I = 1, N
         IF (ABS(X(I)-XMISS).LE.RELXM .OR. ABS(Y(I)-YMISS).LE.RELYM)
     *       GO TO 40
         XI = X(I) - XBAR
         YI = Y(I) - YBAR
         STDX = STDX + XI*XI
         STDY = STDY + YI*YI
         R = R + XI*YI
   40 CONTINUE
      IF (STDX.LE.0.0D0 .OR. STDY.LE.0.0D0) GO TO 220
C
C     AT THIS STAGE XBAR,YBAR,STDX,STDY AND R ARE MEANS,SUMS OF
C     SQUARES AND SUMS OF CROSS-PRODUCTS ABOUT MEANS RESPECTIVELY.
C
      RESULT(21) = FN
      FN = FN - 1.0D0
      RESULT(1) = XBAR
      RESULT(2) = YBAR
      RESULT(3) = SQRT(STDX/FN)
      RESULT(4) = SQRT(STDY/FN)
      RESULT(5) = R/SQRT(STDX*STDY)
      RESULT(13) = 1.0D0
      RESULT(17) = FN - 1.0D0
      RESULT(19) = STDY
      RESULT(20) = FN
      R = R/STDX
      YBAR = YBAR - R*XBAR
      YI = 0.0D0
      DO 60 I = 1, N
         IF (ABS(X(I)-XMISS).LE.RELXM .OR. ABS(Y(I)-YMISS).LE.RELYM)
     *       GO TO 60
         XI = Y(I) - YBAR - R*X(I)
         YI = YI + XI*XI
   60 CONTINUE
C
C     YI IS NOW SUMS OF SQUARES OF DEVIATIONS ABOUT REGRESSION LINE
C
      RESULT(6) = R
      RESULT(7) = YBAR
      RESULT(12) = STDY - YI
      RESULT(14) = RESULT(12)
      RESULT(16) = YI
      YI = YI/RESULT(17)
      RESULT(18) = YI
      IF (YI.LE.1.0D0) GO TO 80
      RESULT(15) = 0.0D0
      IF (ABS(RESULT(14)).GE.YI*SMALL) RESULT(15) = RESULT(14)/YI
      GO TO 100
   80 RESULT(15) = BIGNUM
      IF (ABS(RESULT(14)).LT.YI*BIGNUM) RESULT(15) = RESULT(14)/YI
  100 YBAR = SQRT(YI/STDX)
      XBAR = SQRT(YI*(1.0D0/(FN+1.0D0)+XBAR*XBAR/STDX))
      RESULT(8) = YBAR
      RESULT(9) = XBAR
      IF (YBAR.LE.1.0D0) GO TO 120
      RESULT(10) = 0.0D0
      IF (ABS(RESULT(6)).GE.YBAR*SMALL) RESULT(10) = RESULT(6)/YBAR
      GO TO 140
  120 RESULT(10) = BIGNUM
      IF (ABS(RESULT(6)).LT.YBAR*BIGNUM) RESULT(10) = RESULT(6)/YBAR
  140 IF (XBAR.LE.1.0D0) GO TO 160
      RESULT(11) = 0.0D0
      IF (ABS(RESULT(7)).GE.XBAR*SMALL) RESULT(11) = RESULT(7)/XBAR
      GO TO 180
  160 RESULT(11) = BIGNUM
      IF (ABS(RESULT(7)).LT.XBAR*BIGNUM) RESULT(11) = RESULT(7)/XBAR
C
C     ARRAY RESULT CONTAINS THE FOLLOWING:
C
C     RESULT( 1)     MEAN OF THE INDEPENDENT VARIABLE X
C     RESULT( 2)     MEAN OF THE DEPENDENT VARIABLE Y
C     RESULT( 3)     STANDARD DEVIATION OF THE INDEPENDENT
C     VARIABLE X
C     RESULT( 4)     STANDARD DEVIATION OF THE DEPENDENT
C     VARIABLE Y
C     RESULT( 5)     CORRELATION COEFFICIENT BETWEEN X AND Y
C     RESULT( 6)     REGRESSION COEFFICIENT
C     RESULT( 7)     REGRESSION CONSTANT
C     RESULT( 8)     STANDARD ERROR OF REGRESSION COEFFICIENT
C     RESULT( 9)     STANDARD ERROR OF REGRESSION CONSTANT
C     RESULT(10)     T-STATISTIC FOR REGRESSION COEFFICIENT
C     RESULT(11)     T-STATISTIC FOR REGRESSION CONSTANT
C     RESULT(12)     SUM OF SQUARES DUE TO REGRESSION
C     RESULT(13)     DEGREES OF FREEDOM DUE TO REGRESSION
C     RESULT(14)     MEAN SQUARE DUE TO REGRESSION
C     RESULT(15)     F-STATISTIC
C     RESULT(16)     SUM OF SQUARES ABOUT REGRESSION
C     RESULT(17)     DEGREES OF FREEDOM ABOUT REGRESSION
C     RESULT(18)     MEAN SQUARE ABOUT REGRESSION
C     RESULT(19)     TOTAL SUM OF SQUARES
C     RESULT(20)     TOTAL DEGREES OF FREEDOM
C     RESULT(21)     NUMBER OF CASES USED
C
  180 IFAIL = 0
      RETURN
  200 IERROR = 2
      GO TO 240
  220 IERROR = 3
  240 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END