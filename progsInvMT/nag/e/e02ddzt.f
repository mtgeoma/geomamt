      SUBROUTINE E02DDZ(IOPT,M,X,Y,Z,W,XB,XE,YB,YE,KXX,KYY,S,NXEST,
     *                  NYEST,ETA,TOL,MAXIT,NMAX,KM1,KM2,IB1,IB3,NC,
     *                  INTEST,NREST,NX0,TX,NY0,TY,C,FP,FP0,RANK,FPINT,
     *                  COORD,F,FF,A,Q,BX,BY,SPX,SPY,H,INDEX,ADRES,WRK,
     *                  LWRK,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ETA, FP, FP0, S, TOL, XB, XE, YB, YE
      INTEGER           IB1, IB3, IER, INTEST, IOPT, KM1, KM2, KXX, KYY,
     *                  LWRK, M, MAXIT, NC, NMAX, NREST, NX0, NXEST,
     *                  NY0, NYEST, RANK
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IB1,NC+IB1), BX(NMAX,KM2), BY(NMAX,KM2),
     *                  C(NC), COORD(INTEST), F(NC), FF(NC),
     *                  FPINT(INTEST), H(IB3), Q(IB3,NC+IB3),
     *                  SPX(M,KM1), SPY(M,KM1), TX(NMAX), TY(NMAX),
     *                  W(M), WRK(LWRK), X(M), Y(M), Z(M)
      INTEGER           ADRES(NREST), INDEX(NREST+M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, AJ1, ARG, COS, DMAX, EPS, F1, F2, F3, FAC1,
     *                  FAC2, FPMAX, FPMS, P, P1, P2, P3, PINV, PIV,
     *                  SIGMA, SIN, SQ, STORE, WI, X0, X1, Y0, Y1
      INTEGER           I, I2, IBAND, IBAND1, IBAND3, IBAND4, ICH1,
     *                  ICH3, ICHANG, IFAIL, IN, IROT, ITER, J, JROT,
     *                  JXY, KX, KX1, KX2, KY, KY1, KY2, L, L1, L2, LA,
     *                  LF, LH, LWEST, LX, LY, N, N1, NCOF, NK1X, NK1Y,
     *                  NMINX, NMINY, NREG, NRINT, NUM, NUM1, NX, NXE,
     *                  NXX, NY, NYE, NYY
C     .. Local Arrays ..
      DOUBLE PRECISION  HX(6), HY(6), ZI(1)
C     .. External Functions ..
      DOUBLE PRECISION  E02BEX
      EXTERNAL          E02BEX
C     .. External Subroutines ..
      EXTERNAL          E02BEV, E02BEW, E02DDY, E02ZAF, DROTG, DROT,
     *                  DTBSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, SIGN, SQRT
C     .. Executable Statements ..
C     part 1: Determination of the number of knots and their position.
C     ****************************************************************
C     Given a set of knots we compute the least-squares spline
C     SINF(X,Y), and the corresponding weighted sum of squared
C     residuals FP=F(P=INF). We check whether we can accept the knots:
C     If FP <=S we will continue with the current set of knots.
C     If FP > S we will increase the number of knots and compute the
C      corresponding least-squares spline until finally  FP<=S.
C     The initial choice of knots depends on the value of S and IOPT.
C     If IOPT=0 we first compute the least-squares polynomial of degree
C     KX in X and KY in Y; NX=NMINX=2*KX+2 and NY=NMINY=2*KY+2.
C     FP0=F(0) denotes the corresponding weighted sum of squared
C     residuals
C     If IOPT=1 we start with the knots found at the last call of the
C     routine, except for the case that S>=FP0; then we can compute
C     the least-squares polynomial directly.
C     Eventually the independent variables X and Y (and the
C     corresponding parameters) will be switched if this can reduce
C     the bandwidth of the system to be solved.
C
C     ICHANG denotes whether(1) or not(-1) the directions have been
C     inter-changed.
      ICHANG = -1
      X0 = XB
      X1 = XE
      Y0 = YB
      Y1 = YE
      KX = KXX
      KY = KYY
      KX1 = KX + 1
      KY1 = KY + 1
      NXE = NXEST
      NYE = NYEST
C     Calculation of ACC, the absolute tolerance for the root of F(P)=S.
      ACC = TOL*S
      EPS = SQRT(ETA)
      IF (IOPT.NE.0) THEN
         IF (FP0.GT.S) THEN
            NX = NX0
            NY = NY0
            GO TO 20
         END IF
      END IF
C     Initialization for the least-squares polynomial.
      NMINX = 2*KX1
      NMINY = 2*KY1
      NX = NMINX
      NY = NMINY
      IER = -2
C     Main loop for the different sets of knots. M is a safe upper bound
C     for the number of trials.
   20 DO 780 ITER = 1, M
C        Find the position of the additional knots which are needed
C        for the B-spline representation of S(X,Y).
         DO 40 I = 1, KX1
            TX(I) = X0
            TX(NX-I+1) = X1
   40    CONTINUE
         DO 60 I = 1, KY1
            TY(I) = Y0
            TY(NY-I+1) = Y1
   60    CONTINUE
C        Find NRINT, the total number of knot intervals and NREG,
C        the number of panels in which the approximation domain is
C        subdivided by the intersection of knots.
         NXX = NX - 2*KX1 + 1
         NYY = NY - 2*KY1 + 1
         NRINT = NXX + NYY
         NREG = NXX*NYY
C        Find the bandwidth of the observation matrix A.
C        If necessary, interchange the variables X and Y, in order
C        to obtain a minimal bandwidth.
         IBAND1 = KX*(NY-KY1) + KY
         L = KY*(NX-KX1) + KX
         IF (IBAND1.GT.L) THEN
            IBAND1 = L
            ICHANG = -ICHANG
            DO 80 I = 1, M
               STORE = X(I)
               X(I) = Y(I)
               Y(I) = STORE
   80       CONTINUE
            STORE = X0
            X0 = Y0
            Y0 = STORE
            STORE = X1
            X1 = Y1
            Y1 = STORE
            N = MIN(NX,NY)
            DO 100 I = 1, N
               STORE = TX(I)
               TX(I) = TY(I)
               TY(I) = STORE
  100       CONTINUE
            N1 = N + 1
            IF (NX.LT.NY) THEN
               DO 120 I = N1, NY
                  TX(I) = TY(I)
  120          CONTINUE
            ELSE IF (NX.GT.NY) THEN
               DO 140 I = N1, NX
                  TY(I) = TX(I)
  140          CONTINUE
            END IF
            L = NX
            NX = NY
            NY = L
            L = NXE
            NXE = NYE
            NYE = L
            L = NXX
            NXX = NYY
            NYY = L
            L = KX
            KX = KY
            KY = L
            KX1 = KX + 1
            KY1 = KY + 1
         END IF
         IBAND = IBAND1 + 1
C        Arrange the data points according to the panel they belong to.
         IFAIL = 0
         CALL E02ZAF(NX,NY,TX,TY,M,X,Y,INDEX,NREG+M,ADRES,(NX-7)*(NY-7),
     *               IFAIL)
C        Find NCOF, the number of B-spline coefficients.
         NK1X = NX - KX1
         NK1Y = NY - KY1
         NCOF = NK1X*NK1Y
C        Initialize the observation matrix A.
         DO 160 I = 1, NCOF
            F(I) = ZERO
  160    CONTINUE
         DO 200 I = 1, IBAND
            DO 180 J = IBAND - I + 1, IBAND - I + NCOF
               A(I,J) = ZERO
  180       CONTINUE
  200    CONTINUE
C        Initialize the sum of squared residuals.
         FP = ZERO
C        Fetch the data points in the new order. Main loop for the
C        different panels.
         DO 360 NUM = 1, NREG
C           Fix certain constants for the current panel; JROT records
C           the column number of the first non-zero element in a row
C           of the observation matrix according to a data point of
C           the panel.
            NUM1 = NUM - 1
            LX = NUM1/NYY
            L1 = LX + KX1
            LY = NUM1 - LX*NYY
            L2 = LY + KY1
            JROT = LX*NK1Y + LY
C           Test whether there are still data points in the panel.
            IN = INDEX(NUM+M)
  220       CONTINUE
            IF (IN.NE.0) THEN
C              Fetch a new data point.
               WI = W(IN)
               IF (WI.NE.ZERO) THEN
                  ZI(1) = Z(IN)*WI
C                 Evaluate for the X-direction, the (KX+1) non-zero
C                 B-splines at X(IN).
                  CALL E02BEV(TX,NX,KX,X(IN),L1,HX)
C                 Evaluate for the Y-direction, the (KY+1) non-zero
C                 B-splines at Y(IN).
                  CALL E02BEV(TY,NY,KY,Y(IN),L2,HY)
C                 Store the value of these B-splines in SPX and SPY
C                 respectively.
                  DO 240 I = 1, KX1
                     SPX(IN,I) = HX(I)
  240             CONTINUE
                  DO 260 I = 1, KY1
                     SPY(IN,I) = HY(I)
  260             CONTINUE
C                 Initialize the new row of observation matrix.
                  DO 280 I = 1, IBAND
                     H(I) = ZERO
  280             CONTINUE
C                 Calculate the non-zero elements of the new row by
C                 making the cross products of the non-zero
C                 B-splines in X- and Y-direction.
                  DO 320 I = 1, KX1
                     DO 300 J = 1, KY1
                        H((I-1)*NK1Y+J) = HX(I)*HY(J)*WI
  300                CONTINUE
  320             CONTINUE
C                 Rotate the row into triangle by Givens
C                 transformations .
                  IROT = JROT
                  DO 340 I = 1, IBAND
                     IROT = IROT + 1
                     PIV = ABS(H(I))
                     IF (PIV.NE.ZERO) THEN
C                       Calculate the parameters of the Givens
C                       transformation.
                        AJ1 = ABS(A(IBAND,IROT))
                        CALL DROTG(AJ1,PIV,COS,SIN)
                        COS = SIGN(COS,A(IBAND,IROT))
                        SIN = SIGN(SIN,H(I))
                        A(IBAND,IROT) = AJ1
C                       Apply that transformation to the right hand
C                       side.
                        CALL DROT(1,ZI(1),1,F(IROT),1,COS,-SIN)
C                       Apply that transformation to the left hand side.
                        CALL DROT(IBAND-I,H(I+1),1,A(IBAND-1,IROT+1),
     *                            IB1-1,COS,-SIN)
                     END IF
  340             CONTINUE
C                 Add the contribution of the row to the sum of squares
C                 of residual right hand sides.
                  FP = FP + ZI(1)*ZI(1)
               END IF
C              Find the number of the next data point in the panel.
               IN = INDEX(IN)
               GO TO 220
            END IF
  360    CONTINUE
C        Find DMAX, the maximum value for the diagonal elements in
C        the reduced triangle.
         DMAX = ZERO
         DO 380 I = 1, NCOF
            IF (A(IBAND,I).GT.DMAX) DMAX = A(IBAND,I)
  380    CONTINUE
C        Check whether the observation matrix is rank deficient.
         SIGMA = EPS*DMAX
         DO 400 I = 1, NCOF
            IF (A(IBAND,I).LE.SIGMA) GO TO 460
  400    CONTINUE
C        Backward substitution in case of full rank.
         DO 420 I = 1, NCOF
            C(I) = F(I)
  420    CONTINUE
         CALL DTBSV('U','N','N',NCOF,IBAND-1,A,IB1,C,1)
         RANK = NCOF
         DO 440 I = 1, NCOF
            Q(IBAND,I) = A(IBAND,I)/DMAX
  440    CONTINUE
         GO TO 560
C        In case of rank deficiency, find the minimum norm solution.
C        Check whether there is sufficient working space
  460    LWEST = NCOF*IBAND + NCOF + IBAND
         IF (LWRK.LT.LWEST) THEN
            IER = LWEST
            GO TO 1400
         ELSE
            DO 480 I = 1, NCOF
               FF(I) = F(I)
  480       CONTINUE
            DO 520 I = 1, IBAND
               DO 500 J = IBAND - I + 1, IBAND - I + NCOF
                  Q(I,J) = A(I,J)
  500          CONTINUE
  520       CONTINUE
            LF = 1
            LH = LF + NCOF
            LA = LH + IBAND
            CALL E02DDY(Q,FF,NCOF,IBAND,IB3,SIGMA,C,SQ,RANK,WRK(LA),
     *                  WRK(LF),WRK(LH))
            DO 540 I = 1, NCOF
               Q(IBAND,I) = Q(IBAND,I)/DMAX
  540       CONTINUE
C           Add to the sum of squared residuals, the contribution
C           of reducing the rank.
            FP = FP + SQ
         END IF
  560    IF (IER.EQ.-2) FP0 = FP
C        Test whether the least-squares spline is an acceptable
C        solution.
         FPMS = FP - S
         IF (ABS(FPMS).LE.ACC) THEN
            IF (FP.EQ.ZERO) IER = -1
            IF (NCOF.NE.RANK) IER = -RANK
            GO TO 1400
C           Test whether we can accept the choice of knots.
         ELSE IF (FPMS.LT.ZERO) THEN
            GO TO 800
C           Test whether we cannot further increase the number of knots.
         ELSE IF (NCOF.GT.M) THEN
            IER = 4
            GO TO 1400
         ELSE
            IER = 0
C           Search where to add a new knot.
C           Find for each interval the sum of squared residuals FPINT
C           for the data points having the coordinate belonging to
C           that knot interval. Calculate also COORD which is the same
C           sum, weighted by the position of the data points considered.
            DO 580 I = 1, NRINT
               FPINT(I) = ZERO
               COORD(I) = ZERO
  580       CONTINUE
            DO 660 NUM = 1, NREG
               NUM1 = NUM - 1
               LX = NUM1/NYY
               L1 = LX + 1
               LY = NUM1 - LX*NYY
               L2 = LY + 1 + NXX
               JROT = LX*NK1Y + LY
               IN = INDEX(NUM+M)
  600          CONTINUE
               IF (IN.NE.0) THEN
                  IF (W(IN).NE.ZERO) THEN
                     STORE = ZERO
                     DO 640 I = 1, KX1
                        DO 620 J = 1, KY1
                           STORE = STORE + SPX(IN,I)*SPY(IN,J)*C((I-1)
     *                             *NK1Y+JROT+J)
  620                   CONTINUE
  640                CONTINUE
                     STORE = (W(IN)*(Z(IN)-STORE))**2
                     FPINT(L1) = FPINT(L1) + STORE
                     COORD(L1) = COORD(L1) + STORE*X(IN)
                     FPINT(L2) = FPINT(L2) + STORE
                     COORD(L2) = COORD(L2) + STORE*Y(IN)
                  END IF
                  IN = INDEX(IN)
                  GO TO 600
               END IF
  660       CONTINUE
  680       CONTINUE
C           Find the interval for which FPINT is maximal on the
C           condition that there still can be added a knot.
            L = 0
            FPMAX = ZERO
            L1 = 1
            L2 = NRINT
            IF (NX.EQ.NXE) L1 = NXX + 1
            IF (NY.EQ.NYE) L2 = NXX
            IF (L1.GT.L2) THEN
               IER = 3
               GO TO 1400
            ELSE
               DO 700 I = L1, L2
                  IF (FPMAX.LT.FPINT(I)) THEN
                     L = I
                     FPMAX = FPINT(I)
                  END IF
  700          CONTINUE
C              Test whether we cannot further increase the number
C              of knots.
               IF (L.EQ.0) THEN
                  IER = 5
                  GO TO 1400
               ELSE
C                 Calculate the position of the new knot.
                  ARG = COORD(L)/FPINT(L)
C                 Test in what direction the new knot is going to
C                 be added.
                  IF (L.GT.NXX) THEN
C                    Addition in the Y-direction.
                     JXY = L + KY1 - NXX
                     FPINT(L) = ZERO
                     FAC1 = TY(JXY) - ARG
                     FAC2 = ARG - TY(JXY-1)
                     IF (FAC1.GT.(10*FAC2) .OR. FAC2.GT.(10*FAC1)) THEN
                        GO TO 680
                     ELSE
                        GO TO 740
                     END IF
                  ELSE
C                    Addition in the X-direction.
                     JXY = L + KX1
                     FPINT(L) = ZERO
                     FAC1 = TX(JXY) - ARG
                     FAC2 = ARG - TX(JXY-1)
                     IF (FAC1.GT.(10*FAC2) .OR. FAC2.GT.(10*FAC1))
     *                   GO TO 680
                  END IF
               END IF
            END IF
            DO 720 I = NX, JXY, -1
               TX(I+1) = TX(I)
  720       CONTINUE
            TX(JXY) = ARG
            NX = NX + 1
            GO TO 780
  740       DO 760 I = NY, JXY, -1
               TY(I+1) = TY(I)
  760       CONTINUE
            TY(JXY) = ARG
            NY = NY + 1
         END IF
  780 CONTINUE
C     Restart the computations with the new set of knots.
C     Test whether the least-squares polynomial is a solution of our
C     approximation problem.
  800 IF (IER.EQ.-2) THEN
         GO TO 1400
      ELSE
C        Part 2: Determination of the smoothing spline SP(X,Y)
C        *****************************************************
C        We have determined the number of knots and their position.
C        We now compute the B-spline coefficients of the smoothing
C        spline SP(X,Y). The observation matrix A is extended by
C        the rows of a matrix, expressing that SP(X,Y) must be a
C        polynomial of degree KX in X and KY in Y. The corresponding
C        weights of these additional rows are set to 1./P.
C        Iteratively we than have to determine the value of P
C        such that F(P)=SUM((W(I)*(Z(I)-SP(X(I),Y(I))))**2) be = S.
C        We already know that the least-squares polynomial
C        corresponds to P=0  and that the least-squares spline
C        corresponds to P=infinity. The iteration process which is
C        proposed here makes use of rational interpolation. Since F(P)
C        is a convex and strictly decreasing function of P, it can be
C        approximated by a rational function R(P)= (U*P+V)/(P+W).
C        Three values of P(P1,P2,P3) with corresponding values of
C        F(P) (F1=F(P1)-S,F2=F(P2)-S,F3=F(P3)-S) are used to calculate
C        the new value of P such that R(P)=S. Convergence is
C        guaranteed by taking F1 > 0 and F3 < 0.
C
         KX2 = KX1 + 1
C        Test whether there are interior knots in the X-direction.
C        Evaluate the discontinuity jumps of the KX-th order
C        derivative of the B-splines at the knots TX(L),
C        L=KX+2,...,NX-KX-1.
         IF (NK1X.NE.KX1) CALL E02BEW(TX,NX,KX2,BX,NMAX)
         KY2 = KY1 + 1
C        Test whether there are interior knots in the Y-direction.
C        Evaluate the discontinuity jumps of the KY-th order
C        derivative of the B-splines at the knots TY(L),
C        L=KY+2,...,NY-KY-1.
         IF (NK1Y.NE.KY1) CALL E02BEW(TY,NY,KY2,BY,NMAX)
C        Initial value for P.
         P1 = ZERO
         F1 = FP0 - S
         P3 = -ONE
         F3 = FPMS
         P = ZERO
         DO 820 I = 1, NCOF
            P = P + A(IBAND,I)
  820    CONTINUE
         P = NCOF/P
C        Find the bandwidth of the extended observation matrix.
         IBAND3 = KX1*NK1Y
         IBAND4 = IBAND3 + 1
         ICH1 = 0
         ICH3 = 0
C        Iteration process to find the root of F(P)=S.
         DO 1380 ITER = 1, MAXIT
            PINV = ONE/P
C           Store the triangularized observation matrix into Q.
            DO 840 I = 1, NCOF
               FF(I) = F(I)
  840       CONTINUE
            DO 880 I = 1, IBAND
               DO 860 J = IBAND - I + 1, IBAND - I + NCOF
                  Q(I+IBAND4-IBAND,J) = A(I,J)
  860          CONTINUE
  880       CONTINUE
            DO 920 I = 1, IBAND4 - IBAND
               DO 900 J = IBAND - I + 1, IBAND - I + NCOF
                  Q(I,J) = ZERO
  900          CONTINUE
  920       CONTINUE
            IF (NK1Y.NE.KY1) THEN
C              Extend the observation matrix with the rows of a
C              matrix, expressing that for X=CST, SP(X,Y) must be
C              a polynomial in Y of degree KY.
               DO 1040 I = KY2, NK1Y
                  DO 1020 J = 1, NK1X
C                    Initialize the new row.
                     DO 940 L = 1, IBAND
                        H(L) = ZERO
  940                CONTINUE
C                    Fill in the non-zero elements of the row. JROT
C                    records the column number of the first non-zero
C                    element in the row.
                     DO 960 L = 1, KY2
                        H(L) = BY(I-KY1,L)*PINV
  960                CONTINUE
                     ZI(1) = ZERO
                     JROT = (J-1)*NK1Y + I - KY1
C                    Rotate the new row into triangle by Givens
C                    transformations.
                     DO 1000 IROT = JROT, NCOF
                        PIV = ABS(H(1))
                        I2 = MIN(IBAND1,NCOF-IROT)
                        IF (PIV.NE.ZERO) THEN
C                          Calculate the parameters of the Givens
C                          transformation.
                           AJ1 = ABS(Q(IBAND4,IROT))
                           CALL DROTG(AJ1,PIV,COS,SIN)
                           COS = SIGN(COS,Q(IBAND4,IROT))
                           SIN = SIGN(SIN,H(1))
                           Q(IBAND4,IROT) = AJ1
C                          Apply that Givens transformation to the
C                          right hand side.
                           CALL DROT(1,ZI(1),1,FF(IROT),1,COS,-SIN)
C                          Apply that Givens transformation to the
C                          left hand side.
                           CALL DROT(I2,H(2),1,Q(IBAND4-1,IROT+1),IB3-1,
     *                               COS,-SIN)
                        END IF
                        DO 980 L = 1, I2
                           H(L) = H(L+1)
  980                   CONTINUE
                        H(I2+1) = ZERO
 1000                CONTINUE
 1020             CONTINUE
 1040          CONTINUE
            END IF
            IF (NK1X.NE.KX1) THEN
C              Extend the observation matrix with the rows of a
C              matrix expressing that for Y=CST, SP(X,Y) must be
C              a polynomial in X of degree KX.
               DO 1160 I = KX2, NK1X
                  DO 1140 J = 1, NK1Y
C                    Initialize the new row
                     DO 1060 L = 1, IBAND4
                        H(L) = ZERO
 1060                CONTINUE
C                    Fill in the non-zero elements of the row. JROT
C                    records the column number of the first non-zero
C                    element in the row.
                     DO 1080 L = 1, KX2
                        H((L-1)*NK1Y+1) = BX(I-KX1,L)*PINV
 1080                CONTINUE
                     ZI(1) = ZERO
                     JROT = (I-KX2)*NK1Y + J
C                    Rotate the new row into triangle by Givens
C                    transformations .
                     DO 1120 IROT = JROT, NCOF
                        PIV = ABS(H(1))
                        I2 = MIN(IBAND3,NCOF-IROT)
                        IF (PIV.NE.ZERO) THEN
C                          Calculate the parameters of the Givens
C                          transformation.
                           AJ1 = ABS(Q(IBAND4,IROT))
                           CALL DROTG(AJ1,PIV,COS,SIN)
                           COS = SIGN(COS,Q(IBAND4,IROT))
                           SIN = SIGN(SIN,H(1))
                           Q(IBAND4,IROT) = AJ1
C                          Apply that Givens transformation to the
C                          right hand side.
                           CALL DROT(1,ZI(1),1,FF(IROT),1,COS,-SIN)
C                          Apply that Givens transformation to the
C                          left hand side.
                           CALL DROT(I2,H(2),1,Q(IBAND4-1,IROT+1),IB3-1,
     *                               COS,-SIN)
                        END IF
                        DO 1100 L = 1, I2
                           H(L) = H(L+1)
 1100                   CONTINUE
                        H(I2+1) = ZERO
 1120                CONTINUE
 1140             CONTINUE
 1160          CONTINUE
            END IF
C           Find DMAX, the maximum value for the diagonal elements
C           in the reduced triangle.
            DMAX = ZERO
            DO 1180 I = 1, NCOF
               IF (Q(IBAND4,I).GT.DMAX) DMAX = Q(IBAND4,I)
 1180       CONTINUE
C           Check whether the matrix is rank deficient.
            SIGMA = EPS*DMAX
            DO 1200 I = 1, NCOF
               IF (Q(IBAND4,I).LE.SIGMA) GO TO 1240
 1200       CONTINUE
C           Backward substitution in case of full rank.
            DO 1220 I = 1, NCOF
               C(I) = FF(I)
 1220       CONTINUE
            CALL DTBSV('U','N','N',NCOF,IBAND4-1,Q,IB3,C,1)
            RANK = NCOF
            GO TO 1260
C           In case of rank deficiency, find the minimum norm solution.
 1240       LWEST = NCOF*IBAND4 + NCOF + IBAND4
            IF (LWRK.LT.LWEST) THEN
               IER = LWEST
               GO TO 1400
            ELSE
               LF = 1
               LH = LF + NCOF
               LA = LH + IBAND4
               CALL E02DDY(Q,FF,NCOF,IBAND4,IB3,SIGMA,C,SQ,RANK,WRK(LA),
     *                     WRK(LF),WRK(LH))
            END IF
 1260       DO 1280 I = 1, NCOF
               Q(IBAND4,I) = Q(IBAND4,I)/DMAX
 1280       CONTINUE
C           Compute F(P).
            FP = ZERO
            DO 1360 NUM = 1, NREG
               NUM1 = NUM - 1
               LX = NUM1/NYY
               LY = NUM1 - LX*NYY
               JROT = LX*NK1Y + LY
               IN = INDEX(NUM+M)
 1300          CONTINUE
               IF (IN.NE.0) THEN
                  IF (W(IN).NE.ZERO) THEN
                     STORE = ZERO
                     DO 1340 I = 1, KX1
                        DO 1320 J = 1, KY1
                           STORE = STORE + SPX(IN,I)*SPY(IN,J)*C((I-1)
     *                             *NK1Y+JROT+J)
 1320                   CONTINUE
 1340                CONTINUE
                     FP = FP + (W(IN)*(Z(IN)-STORE))**2
                  END IF
                  IN = INDEX(IN)
                  GO TO 1300
               END IF
 1360       CONTINUE
C           Test whether the approximation SP(X,Y) is an acceptable
C           solution.
            FPMS = FP - S
            IF (ABS(FPMS).LE.ACC) THEN
               IF (NCOF.NE.RANK) IER = -RANK
               GO TO 1400
C              Test whether the maximum allowable number of iterations
C              has been reached.
            ELSE IF (ITER.EQ.MAXIT) THEN
               IER = 6
               GO TO 1400
            ELSE
C              Carry out one more step of the iteration process.
               P2 = P
               F2 = FPMS
               IF (ICH3.EQ.0) THEN
                  IF ((F2-F3).LE.ACC) THEN
C                    Our initial choice of P is too large.
                     P3 = P2
                     F3 = F2
                     P = P/25
                     IF (P.LE.P1) P = (P1*9+P2)/10
                     GO TO 1380
                  ELSE IF (F2.LT.ZERO) THEN
                     ICH3 = 1
                  END IF
               END IF
               IF (ICH1.EQ.0) THEN
                  IF ((F1-F2).LE.ACC) THEN
C                    Our initial choice of P is too small
                     P1 = P2
                     F1 = F2
                     P = P*25
                     IF (P3.GE.ZERO) THEN
                        IF (P.GE.P3) P = (P2+P3*9)/10
                     END IF
                     GO TO 1380
                  ELSE IF (F2.GT.ZERO) THEN
                     ICH1 = 1
                  END IF
               END IF
C              Test whether the iteration process proceeds as
C              theoretically expected.
               IF (F2.GE.F1 .OR. F2.LE.F3) THEN
                  IER = 6
                  GO TO 1400
               ELSE
C                 Find the new value of P.
                  P = E02BEX(P1,F1,P2,F2,P3,F3)
               END IF
            END IF
 1380    CONTINUE
      END IF
C     Test whether X and Y are in the original order.
 1400 IF (ICHANG.EQ.1) THEN
C        Interchange X and Y once more.
         DO 1440 I = 1, NK1X
            DO 1420 J = 1, NK1Y
               F((J-1)*NK1X+I) = C((I-1)*NK1Y+J)
 1420       CONTINUE
 1440    CONTINUE
         DO 1460 I = 1, NCOF
            C(I) = F(I)
 1460    CONTINUE
         DO 1480 I = 1, M
            STORE = X(I)
            X(I) = Y(I)
            Y(I) = STORE
 1480    CONTINUE
         N = MIN(NX,NY)
         DO 1500 I = 1, N
            STORE = TX(I)
            TX(I) = TY(I)
            TY(I) = STORE
 1500    CONTINUE
         N1 = N + 1
         IF (NX.LT.NY) THEN
            DO 1520 I = N1, NY
               TX(I) = TY(I)
 1520       CONTINUE
         ELSE IF (NX.GT.NY) THEN
            DO 1540 I = N1, NX
               TY(I) = TX(I)
 1540       CONTINUE
         END IF
         L = NX
         NX = NY
         NY = L
      END IF
      NX0 = NX
      NY0 = NY
C     Ignore negative IER for NAG version.
      IF (IER.LT.0) IER = 0
      RETURN
      END
