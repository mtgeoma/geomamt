      SUBROUTINE G01DDF(X,N,CALWTS,A,W,PW,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     MODIFICATION OF ALGORITHM AS 181  APPL. STATIST. (1982)
C
C     CALCULATES SHAPIRO AND WILK W STATISTIC AND ITS SIG. LEVEL
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01DDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PW, W
      INTEGER           IFAIL, N
      LOGICAL           CALWTS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AL, AN, EPS, EU3, FIVE, LAMDA, ONE, ONEPT4, PI6,
     *                  SDY, SSQ, STQR, T2, T3, THREE, TQR, UN, WW, X1,
     *                  X2, XM, Y, YBAR, Z, ZERO
      INTEGER           I, IERR, J, N2, N3, NC
      LOGICAL           UP
C     .. Local Arrays ..
      DOUBLE PRECISION  C(0:5), C1(5,3), C2(5,3), UNH(3), UNL(3),
     *                  WA(0:2), WB(0:3), WC(0:3), WD(0:5), WE(0:5),
     *                  WF(0:6)
      INTEGER           NC1(3), NC2(3)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01DDZ, S15ACF
      INTEGER           P01ABF
      EXTERNAL          G01DDZ, S15ACF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01DDY
C     .. Intrinsic Functions ..
      INTRINSIC         ASIN, EXP, LOG, SQRT, DBLE
C     .. Data statements ..
      DATA              WA(0), WA(1), WA(2)/0.118898D0, 0.133414D0,
     *                  0.327907D0/, WB(0), WB(1), WB(2),
     *                  WB(3)/-0.37542D0, -0.492145D0, -1.124332D0,
     *                  -0.199422D0/, WC(0), WC(1), WC(2),
     *                  WC(3)/-3.15805D0, 0.729399D0, 3.01855D0,
     *                  1.558776D0/, WD(0), WD(1), WD(2), WD(3), WD(4),
     *                  WD(5)/0.480385D0, 0.318828D0, 0.0D0,
     *                  -0.0241665D0, 0.00879701D0, 0.002989646D0/,
     *                  WE(0), WE(1), WE(2), WE(3), WE(4),
     *                  WE(5)/-1.91487D0, -1.37888D0, -0.04183209D0,
     *                  0.1066339D0, -0.03513666D0, -0.01504614D0/,
     *                  WF(0), WF(1), WF(2), WF(3), WF(4), WF(5),
     *                  WF(6)/-3.73538D0, -1.015807D0, -0.331885D0,
     *                  0.1773538D0, -0.01638782D0, -0.03215018D0,
     *                  0.003852646D0/
      DATA              C1(1,1), C1(2,1), C1(3,1), C1(4,1), C1(5,1),
     *                  C1(1,2), C1(2,2), C1(3,2), C1(4,2), C1(5,2),
     *                  C1(1,3), C1(2,3), C1(3,3), C1(4,3),
     *                  C1(5,3)/-1.26233D0, 1.87969D0, 0.0649583D0,
     *                  -0.0475604D0, -0.0139682D0, -2.28135D0,
     *                  2.26186D0, 0.0D0, 0.0D0, -0.00865763D0,
     *                  -3.30623D0, 2.76287D0, -0.83484D0, 1.20857D0,
     *                  -0.507590D0/
      DATA              C2(1,1), C2(2,1), C2(3,1), C2(4,1), C2(5,1),
     *                  C2(1,2), C2(2,2), C2(3,2), C2(4,2), C2(5,2),
     *                  C2(1,3), C2(2,3), C2(3,3), C2(4,3),
     *                  C2(5,3)/-0.287696D0, 1.78953D0, -0.180114D0,
     *                  0.0D0, 0.0D0, -1.63638D0, 5.60924D0, -3.63738D0,
     *                  1.08439D0, 0.0D0, -5.991908D0, 21.04575D0,
     *                  -24.58061D0, 13.78661D0, -2.835295D0/
      DATA              UNL(1), UNL(2), UNL(3)/-3.8D0, -3.0D0, -1.0D0/,
     *                  UNH(1), UNH(2), UNH(3)/8.6D0, 5.8D0, 5.4D0/
      DATA              NC1(1), NC1(2), NC1(3)/5, 5, 5/, NC2(1), NC2(2),
     *                  NC2(3)/3, 4, 5/
      DATA              PI6/1.90985932D0/, STQR/1.04719755D0/,
     *                  ZERO/0.0D0/, TQR/0.75D0/, ONE/1.0D0/,
     *                  ONEPT4/1.4D0/, THREE/3.0D0/, FIVE/5.0D0/
C     .. Executable Statements ..
      IF (N.LE.2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) '** ON ENTRY, THE VALUE OF N (', N,
     *     ' ) IS LESS THAN 3'
         GO TO 280
      END IF
      IF (N.GT.2000) THEN
         IERR = 2
         WRITE (P01REC,FMT=99999) '** ON ENTRY, THE VALUE OF N (', N,
     *     ' ) IS GREATER THAN 2000'
         GO TO 280
      END IF
C
C     CHECK WHETHER ALL OBSERVATIONS ARE EQUAL
C
      IERR = 3
      IF (X(1).EQ.X(N)) THEN
         WRITE (P01REC,FMT=99999)
     *     '** ON ENTRY, ALL THE ELEMENTS IN X ARE EQUAL'
         GO TO 280
      END IF
      IF (X(1).LT.X(N)) THEN
         UP = .TRUE.
      ELSE
         UP = .FALSE.
      END IF
C
C     CHECK WHETHER DATA IS SORTED INTO EITHER ASCENDING OR
C     DESCENDING ORDER AND CALCULATE THE SUM OF SQUARES ABOUT
C     THE MEAN USING WEST'S ALGORITHM
C
      XM = X(1)
      SSQ = 0.0D0
      DO 40 I = 2, N
C
         X1 = X(I-1)
         X2 = X(I)
         IF (UP) GO TO 20
         X1 = X(I)
         X2 = X(I-1)
   20    IF (X2.LT.X1) THEN
            WRITE (P01REC,FMT=99998)
     *        '** ON ENTRY, THE ELEMENTS OF X ARE NOT IN ASCENDING OR ',
     *        'DESCENDING ORDER'
            GO TO 280
         END IF
         T2 = X(I) - XM
         T3 = T2/DBLE(I)
         XM = XM + T3
         SSQ = SSQ + DBLE(I-1)*T2*T3
C
   40 CONTINUE
C
C     NO ERRORS IN INPUT DATA SO SET IFAIL = 0
C
      IFAIL = 0
      IF (CALWTS) THEN
         CALL G01DDY(N,A,EPS)
      ELSE
         EPS = A(1)*A(1)/(ONE-ONE/DBLE(N))
      END IF
C
C     CALCULATE W
C
      N2 = N/2
      W = ZERO
      PW = ONE
      AN = N
      I = N
      DO 60 J = 1, N2
         W = W + A(J)*(X(I)-X(J))
         I = I - 1
   60 CONTINUE
      W = W*W/SSQ
      IF (W.LT.ONE) GO TO 80
C
C     WHEN W IS .GE. 1.0 RETURN W AS 1.0 AND THE SIGNIFICANCE
C     LEVEL PW AS 1.0
C
      W = ONE
      RETURN
C
C     GET SIGNIFICANCE LEVEL OF W
C
   80 IF (N.LE.6) GO TO 140
C
C     N BETWEEN 7 AND 2000 ... TRANSFORM W TO Y, GET MEAN AND SD,
C     STANDARDIZE AND GET SIGNIFICANCE LEVEL
C
      IF (N.GT.20) GO TO 100
      AL = LOG(AN) - THREE
      LAMDA = G01DDZ(WA,2,AL)
      YBAR = EXP(G01DDZ(WB,3,AL))
      SDY = EXP(G01DDZ(WC,3,AL))
      GO TO 120
  100 AL = LOG(AN) - FIVE
      LAMDA = G01DDZ(WD,5,AL)
      YBAR = EXP(G01DDZ(WE,5,AL))
      SDY = EXP(G01DDZ(WF,6,AL))
  120 Y = (ONE-W)**LAMDA
      Z = (Y-YBAR)/SDY
      PW = S15ACF(Z,IFAIL)
      RETURN
C
C     DEAL WITH N LESS THAN 7 (EXACT SIGNIFICANCE LEVEL FOR N=3).
C
  140 IF (W.LE.EPS) GO TO 260
      WW = W
      IF (N.EQ.3) GO TO 240
      UN = LOG((W-EPS)/(ONE-W))
      N3 = N - 3
      IF (UN.LT.UNL(N3)) GO TO 260
      IF (UN.GE.ONEPT4) GO TO 180
      NC = NC1(N3)
      DO 160 I = 1, NC
         C(I-1) = C1(I,N3)
  160 CONTINUE
      EU3 = EXP(G01DDZ(C,NC-1,UN))
      GO TO 220
C
C     IF UN IS GREATER THAN ITS UPPER BOUND RETURN PW AS 1.0
C
  180 IF (UN.GT.UNH(N3)) RETURN
      NC = NC2(N3)
      DO 200 I = 1, NC
         C(I-1) = C2(I,N3)
  200 CONTINUE
      UN = LOG(UN)
      EU3 = EXP(EXP(G01DDZ(C,NC-1,UN)))
  220 WW = (EU3+TQR)/(ONE+EU3)
  240 PW = PI6*(ASIN(SQRT(WW))-STQR)
      RETURN
C
C     IF W IS LESS THAN OR EQUAL TO ITS SMALLEST POSSIBLE VALUE
C     OR UN IS LESS THAN ITS LOWER BOUND RETURN PW AS 0.0
C
  260 PW = ZERO
      RETURN
  280 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,P01REC)
      RETURN
C
99999 FORMAT (1X,A,I16,A)
99998 FORMAT (1X,A,A)
      END
