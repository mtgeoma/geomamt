      SUBROUTINE D01APX(F,A,B,BL,BR,ALFA,BETA,RI,RJ,RG,RH,RESULT,ABSERR,
     *                  RESASC,INTEGR,NEV)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QC25S.
C     ..................................................................
C
C        PURPOSE
C           TO COMPUTE I = INTEGRAL OF F*W OVER (BL,BR), WITH ERROR
C           ESTIMATE, WHERE THE WEIGHT FUNCTION W HAS A SINGULAR
C           BEHAVIOUR OF ALGEBRAICO-LOGARITHMIC TYPE AT THE POINTS
C           A AND/OR B. (BL,BR) IS A PART OF (A,B).
C
C        PARAMETERS
C           F      - REAL
C                    FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
C                    F(X). THE ACTUAL NAME FOR F NEEDS TO BE DECLARED
C                    E X T E R N A L  IN THE DRIVER PROGRAM.
C
C           A      - REAL
C                    LEFT END POINT OF THE ORIGINAL INTERVAL
C
C           B      - REAL
C                    RIGHT END POINT OF THE ORIGINAL INTERVAL, B.GT.A
C
C           BL     - REAL
C                    LOWER LIMIT OF INTEGRATION, BL.GE.A
C
C           BR     - REAL
C                    UPPER LIMIT OF INTEGRATION, BR.LE.B
C
C           ALFA   - REAL
C                    PARAMETER IN THE WEIGHT FUNCTION
C
C           BETA   - REAL
C                    PARAMETER IN THE WEIGHT FUNCTION
C
C           RI,RJ,RG,RH - REAL
C                    MODIFIED CHEBYSHEV MOMENTS FOR THE APPLICATION
C                    OF THE GENERALIZED CLENSHAW-CURTIS METHOD
C                    (COMPUTED IN SUBROUTINE D01APW)
C
C           RESULT - REAL
C                    APPROXIMATION TO THE INTEGRAL
C                    RESULT IS COMPUTED BY USING A GENERALIZED
C                    CLENSHAW-CURTIS METHOD IF B1 = A OR BR = B.
C                    IN ALL OTHER CASES THE 15-POINT KRONROD RULE IS
C                    APPLIED, OBTAINED BY OPTIMAL ADDITION OF ABSCISSAE
C                    TO THE 7-POINT GAUSS RULE.
C
C           ABSERR - REAL
C                    ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
C                    WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)
C
C           RESASC - REAL
C                    APPROXIMATION TO THE INTEGRAL OF ABS(F*W-I/(B-A))
C
C           INTEGR - INTEGER
C                    WHICH DETERMINES THE WEIGHT FUNCTION
C                    = 1   W(X) = (X-A)**ALFA*(B-X)**BETA
C                    = 2   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
C                    = 3   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
C                    = 4   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*
C                                 LOG(B-X)
C
C           NEV    - INTEGER
C                    NUMBER OF INTEGRAND EVALUATIONS
C
C     ..................................................................
C
C           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
C           K = 1, ..., 11, TO BE USED FOR THE COMPUTATION OF THE
C           CHEBYSHEV SERIES EXPANSION OF F.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, ALFA, B, BETA, BL, BR, RESASC, RESULT
      INTEGER           INTEGR, NEV
C     .. Array Arguments ..
      DOUBLE PRECISION  RG(25), RH(25), RI(25), RJ(25)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  CENTR, DC, FACTOR, FIX, HLGTH, RES12, RES24,
     *                  RESABS, U
      INTEGER           I, ISYM
C     .. Local Arrays ..
      DOUBLE PRECISION  CHEB12(13), CHEB24(25), FVAL(25), X(11)
C     .. External Functions ..
      DOUBLE PRECISION  D01APY
      EXTERNAL          D01APY
C     .. External Subroutines ..
      EXTERNAL          D01ANX, D01APZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG
C     .. Data statements ..
      DATA              X(1)/0.991444861373810411144557526928563D+00/,
     *                  X(2)/0.965925826289068286749743199728897D+00/,
     *                  X(3)/0.923879532511286756128183189396788D+00/,
     *                  X(4)/0.866025403784438646763723170752936D+00/,
     *                  X(5)/0.793353340291235164579776961501299D+00/,
     *                  X(6)/0.707106781186547524400844362104849D+00/,
     *                  X(7)/0.608761429008720639416097542898164D+00/,
     *                  X(8)/0.500000000000000000000000000000000D+00/,
     *                  X(9)/0.382683432365089771728459984030399D+00/,
     *                  X(10)/0.258819045102520762348898837624048D+00/,
     *                  X(11)/0.130526192220051591548406227895489D+00/
C     .. Executable Statements ..
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
C                    (BR-BL)*0.5*COS(K*PI/24)+(BR+BL)*0.5
C                    K = 0, ..., 24
C           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
C                    OF DEGREE 12, FOR THE FUNCTION F, IN THE INTERVAL
C                    (BL,BR)
C           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
C                    OF DEGREE 24, FOR THE FUNCTION F, IN THE INTERVAL
C                    (BL,BR)
C           RES12  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB12
C           RES24  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB24
C           D01APY  - EXTERNAL FUNCTION SUBPROGRAM DEFINING THE FOUR
C                    POSSIBLE WEIGHT FUNCTIONS
C           HLGTH  - HALF-LENGTH OF THE INTERVAL (BL,BR)
C           CENTR  - MID POINT OF THE INTERVAL (BL,BR)
C
      NEV = 25
      IF (BL.EQ.A .AND. (ALFA.NE.0.0D+00 .OR. INTEGR.EQ.2 .OR.
     *    INTEGR.EQ.4)) GO TO 20
      IF (BR.EQ.B .AND. (BETA.NE.0.0D+00 .OR. INTEGR.EQ.3 .OR.
     *    INTEGR.EQ.4)) GO TO 280
C
C           IF A.GT.BL AND B.LT.BR, APPLY THE 15-POINT GAUSS-KRONROD
C           SCHEME.
C
      CALL D01APZ(F,D01APY,A,B,ALFA,BETA,INTEGR,BL,BR,RESULT,ABSERR,
     *            RESABS,RESASC)
      NEV = 15
      GO TO 540
C
C           THIS PART OF THE PROGRAM IS EXECUTED ONLY IF A = BL.
C           ----------------------------------------------------
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F1 = (0.5*(B+B-BR-A)-0.5*(BR-A)*X)**BETA
C                  *F(0.5*(BR-A)*X+0.5*(BR+A))
C
   20 HLGTH = 5.0D-01*(BR-BL)
      CENTR = 5.0D-01*(BR+BL)
      FIX = B - CENTR
      FVAL(1) = 5.0D-01*F(HLGTH+CENTR)*(FIX-HLGTH)**BETA
      FVAL(13) = F(CENTR)*(FIX**BETA)
      FVAL(25) = 5.0D-01*F(CENTR-HLGTH)*(FIX+HLGTH)**BETA
      DO 40 I = 2, 12
         U = HLGTH*X(I-1)
         ISYM = 26 - I
         FVAL(I) = F(U+CENTR)*(FIX-U)**BETA
         FVAL(ISYM) = F(CENTR-U)*(FIX+U)**BETA
   40 CONTINUE
      FACTOR = HLGTH**(ALFA+1.0D+00)
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      RES12 = 0.0D+00
      RES24 = 0.0D+00
      IF (INTEGR.GT.2) GO TO 140
      CALL D01ANX(X,FVAL,CHEB12,CHEB24)
C
C           INTEGR = 1  (OR 2)
C
      DO 60 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RI(I)
         RES24 = RES24 + CHEB24(I)*RI(I)
   60 CONTINUE
      DO 80 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RI(I)
   80 CONTINUE
      IF (INTEGR.EQ.1) GO TO 260
C
C           INTEGR = 2
C
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0D+00
      RES24 = 0.0D+00
      DO 100 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RG(I)
         RES24 = RES24 + CHEB24(I)*RG(I)
  100 CONTINUE
      DO 120 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RG(I)
  120 CONTINUE
      GO TO 260
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F4 = F1*LOG(0.5*(B+B-BR-A)-0.5*(BR-A)*X)
C
  140 FVAL(1) = FVAL(1)*LOG(FIX-HLGTH)
      FVAL(13) = FVAL(13)*LOG(FIX)
      FVAL(25) = FVAL(25)*LOG(FIX+HLGTH)
      DO 160 I = 2, 12
         U = HLGTH*X(I-1)
         ISYM = 26 - I
         FVAL(I) = FVAL(I)*LOG(FIX-U)
         FVAL(ISYM) = FVAL(ISYM)*LOG(FIX+U)
  160 CONTINUE
      CALL D01ANX(X,FVAL,CHEB12,CHEB24)
C
C           INTEGR = 3  (OR 4)
C
      DO 180 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RI(I)
         RES24 = RES24 + CHEB24(I)*RI(I)
  180 CONTINUE
      DO 200 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RI(I)
  200 CONTINUE
      IF (INTEGR.EQ.3) GO TO 260
C
C           INTEGR = 4
C
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0D+00
      RES24 = 0.0D+00
      DO 220 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RG(I)
         RES24 = RES24 + CHEB24(I)*RG(I)
  220 CONTINUE
      DO 240 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RG(I)
  240 CONTINUE
  260 RESULT = (RESULT+RES24)*FACTOR
      ABSERR = (ABSERR+ABS(RES24-RES12))*FACTOR
      GO TO 540
C
C           THIS PART OF THE PROGRAM IS EXECUTED ONLY IF B = BR.
C           ----------------------------------------------------
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F2 = (0.5*(B+BL-A-A)+0.5*(B-BL)*X)**ALFA
C                *F(0.5*(B-BL)*X+0.5*(B+BL))
C
  280 HLGTH = 5.0D-01*(BR-BL)
      CENTR = 5.0D-01*(BR+BL)
      FIX = CENTR - A
      FVAL(1) = 5.0D-01*F(HLGTH+CENTR)*(FIX+HLGTH)**ALFA
      FVAL(13) = F(CENTR)*(FIX**ALFA)
      FVAL(25) = 5.0D-01*F(CENTR-HLGTH)*(FIX-HLGTH)**ALFA
      DO 300 I = 2, 12
         U = HLGTH*X(I-1)
         ISYM = 26 - I
         FVAL(I) = F(U+CENTR)*(FIX+U)**ALFA
         FVAL(ISYM) = F(CENTR-U)*(FIX-U)**ALFA
  300 CONTINUE
      FACTOR = HLGTH**(BETA+1.0D+00)
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      RES12 = 0.0D+00
      RES24 = 0.0D+00
      IF (INTEGR.EQ.2 .OR. INTEGR.EQ.4) GO TO 400
C
C           INTEGR = 1  (OR 3)
C
      CALL D01ANX(X,FVAL,CHEB12,CHEB24)
      DO 320 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RJ(I)
         RES24 = RES24 + CHEB24(I)*RJ(I)
  320 CONTINUE
      DO 340 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RJ(I)
  340 CONTINUE
      IF (INTEGR.EQ.1) GO TO 520
C
C           INTEGR = 3
C
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0D+00
      RES24 = 0.0D+00
      DO 360 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RH(I)
         RES24 = RES24 + CHEB24(I)*RH(I)
  360 CONTINUE
      DO 380 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RH(I)
  380 CONTINUE
      GO TO 520
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F3 = F2*LOG(0.5*(B-BL)*X+0.5*(B+BL-A-A))
C
  400 FVAL(1) = FVAL(1)*LOG(HLGTH+FIX)
      FVAL(13) = FVAL(13)*LOG(FIX)
      FVAL(25) = FVAL(25)*LOG(FIX-HLGTH)
      DO 420 I = 2, 12
         U = HLGTH*X(I-1)
         ISYM = 26 - I
         FVAL(I) = FVAL(I)*LOG(U+FIX)
         FVAL(ISYM) = FVAL(ISYM)*LOG(FIX-U)
  420 CONTINUE
      CALL D01ANX(X,FVAL,CHEB12,CHEB24)
C
C           INTEGR = 2  (OR 4)
C
      DO 440 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RJ(I)
         RES24 = RES24 + CHEB24(I)*RJ(I)
  440 CONTINUE
      DO 460 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RJ(I)
  460 CONTINUE
      IF (INTEGR.EQ.2) GO TO 520
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0D+00
      RES24 = 0.0D+00
C
C           INTEGR = 4
C
      DO 480 I = 1, 13
         RES12 = RES12 + CHEB12(I)*RH(I)
         RES24 = RES24 + CHEB24(I)*RH(I)
  480 CONTINUE
      DO 500 I = 14, 25
         RES24 = RES24 + CHEB24(I)*RH(I)
  500 CONTINUE
  520 RESULT = (RESULT+RES24)*FACTOR
      ABSERR = (ABSERR+ABS(RES24-RES12))*FACTOR
  540 RETURN
      END
