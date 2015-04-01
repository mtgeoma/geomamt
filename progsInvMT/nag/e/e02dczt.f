      SUBROUTINE E02DCZ(IOPT,X,MX,Y,MY,Z,XB,XE,YB,YE,KX,KY,S,NXEST,
     *                  NYEST,TOL,MAXIT,NC,NX,TX,NY,TY,C,FP,FP0,FPOLD,
     *                  REDUCX,REDUCY,FPINTX,FPINTY,LASTDI,NPLUSX,
     *                  NPLUSY,NRX,NRY,NRDATX,NRDATY,WRK,LWRK,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16A REVISED. IER-983 (JUN 1993).
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FP, FP0, FPOLD, REDUCX, REDUCY, S, TOL, XB, XE,
     *                  YB, YE
      INTEGER           IER, IOPT, KX, KY, LASTDI, LWRK, MAXIT, MX, MY,
     *                  NC, NPLUSX, NPLUSY, NX, NXEST, NY, NYEST
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NC), FPINTX(NXEST), FPINTY(NYEST), TX(NXEST),
     *                  TY(NYEST), WRK(LWRK), X(MX), Y(MY), Z(MX*MY)
      INTEGER           NRDATX(NXEST), NRDATY(NYEST), NRX(MX), NRY(MY)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, F1, F2, F3, FPMS, P, P1, P2, P3
      INTEGER           I, ICH1, ICH3, IFBX, IFBY, IFSX, IFSY, ITER, J,
     *                  K3, KX1, KX2, KY1, KY2, L, LAX, LAY, LBX, LBY,
     *                  LQ, LRI, LSX, LSY, MK1, MM, MPM, MYNX, NCOF,
     *                  NK1X, NK1Y, NMAXX, NMAXY, NMINX, NMINY, NPL1,
     *                  NPLX, NPLY, NRINTX, NRINTY, NXE, NXK, NYE
C     .. External Functions ..
      DOUBLE PRECISION  E02BEX
      EXTERNAL          E02BEX
C     .. External Subroutines ..
      EXTERNAL          E02BEY, E02DCY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C     We partition the working space.
      KX1 = KX + 1
      KY1 = KY + 1
      KX2 = KX1 + 1
      KY2 = KY1 + 1
      LSX = 1
      LSY = LSX + MX*KX1
      LRI = LSY + MY*KY1
      MM = MAX(NXEST,MY)
      LQ = LRI + MM
      MYNX = NXEST*MY
      LAX = LQ + MYNX
      NXK = NXEST*KX2
      LBX = LAX + NXK + KX2*KX2
      LAY = LBX + NXK
      LBY = LAY + (NYEST+KY2)*KY2
C     Part 1: Determination of the number of knots and their position.
C     ****************************************************************
C     Given a set of knots we compute the least-squares spline
C     SINF(X,Y), and the corresponding sum of squared residuals
C     FP=F(P=INF). We check whether we can accept the knots:
C     if FP <=S we will continue with the current set of knots.
C     if FP > S we will increase the number of knots and compute the
C       corresponding least-squares spline until finally FP<=S.
C     The initial choice of knots depends on the value of S and IOPT.
C     If S=0 we have spline interpolation; in that case the number of
C     knots equals NMAXX = MX+KX+1  and  NMAXY = MY+KY+1.
C     If S>0 and
C     *IOPT=0 we first compute the least-squares polynomial of degree
C      KX in X and KY in Y; NX=NMINX=2*KX+2 and NY=NYMIN=2*KY+2.
C     *IOPT=1 we start with the knots found at the last call of the
C      routine, except for the case that S > FP0; then we can compute
C      the least-squares polynomial directly.
C
C     Determine the number of knots for polynomial approximation.
      NMINX = 2*KX1
      NMINY = 2*KY1
C     ACC denotes the absolute tolerance for the root of F(P)=S.
      ACC = TOL*S
C     Find NMAXX and NMAXY which denote the number of knots in X- and Y-
C     direction in case of spline interpolation.
      NMAXX = MX + KX1
      NMAXY = MY + KY1
C     Find NXE and NYE which denote the maximum number of knots
C     allowed in each direction
      NXE = MIN(NMAXX,NXEST)
      NYE = MIN(NMAXY,NYEST)
      IF (S.GT.ZERO) THEN
C        If S > 0 our initial choice of knots depends on the value
C        of IOPT.
         IF (IOPT.NE.0) THEN
            IF (FP0.GT.S) THEN
C              If IOPT=1 and FP0 > S we start computing the
C              least- squares spline according to the set of knots
C              found at the last call of the routine. We determine
C              the number of grid coordinates X(I) inside each knot
C              interval (TX(L),TX(L+1)).
               L = KX2
               J = 1
               NRDATX(1) = 0
               MPM = MX - 1
               DO 20 I = 2, MPM
                  NRDATX(J) = NRDATX(J) + 1
                  IF (X(I).GE.TX(L)) THEN
                     NRDATX(J) = NRDATX(J) - 1
                     L = L + 1
                     J = J + 1
                     NRDATX(J) = 0
                  END IF
   20          CONTINUE
C              We determine the number of grid coordinates Y(I) inside
C              each knot interval (TY(L),TY(L+1)).
               L = KY2
               J = 1
               NRDATY(1) = 0
               MPM = MY - 1
               DO 40 I = 2, MPM
                  NRDATY(J) = NRDATY(J) + 1
                  IF (Y(I).GE.TY(L)) THEN
                     NRDATY(J) = NRDATY(J) - 1
                     L = L + 1
                     J = J + 1
                     NRDATY(J) = 0
                  END IF
   40          CONTINUE
               GO TO 140
            END IF
         END IF
C        If IOPT=0 or IOPT=1 and S>=FP0, we start computing the
C        least-squares polynomial of degree KX in X and KY in Y
C        (which is a spline without interior knots).
         NX = NMINX
         NY = NMINY
         NRDATX(1) = MX - 2
         NRDATY(1) = MY - 2
         LASTDI = 0
         NPLUSX = 0
         NPLUSY = 0
         FP0 = ZERO
         FPOLD = ZERO
         REDUCX = ZERO
         REDUCY = ZERO
      ELSE
C        If S = 0, S(X,Y) is an interpolating spline.
         NX = NMAXX
         NY = NMAXY
C        Test whether the required storage space exceeds the
C        available one.
         IF (NY.GT.NYEST .OR. NX.GT.NXEST) THEN
            IER = 4
            GO TO 360
         ELSE
C           Find the position of the interior knots in case of
C           interpolation. The knots in the X-direction.
            MK1 = MX - KX1
            IF (MK1.NE.0) THEN
               K3 = KX/2
               IF (K3*2.EQ.KX) THEN
                  DO 60 L = 1, MK1
                     TX(KX1+L) = (X(K3+L+1)-X(K3+L))/2
   60             CONTINUE
               ELSE
                  DO 80 L = 1, MK1
                     TX(KX1+L) = X(K3+L+1)
   80             CONTINUE
               END IF
            END IF
C           The knots in the Y-direction.
            MK1 = MY - KY1
            IF (MK1.NE.0) THEN
               K3 = KY/2
               IF (K3*2.EQ.KY) THEN
                  DO 100 L = 1, MK1
                     TY(KY1+L) = (Y(K3+L+1)-Y(K3+L))/2
  100             CONTINUE
               ELSE
                  DO 120 L = 1, MK1
                     TY(KY1+L) = Y(K3+L+1)
  120             CONTINUE
               END IF
            END IF
         END IF
      END IF
  140 MPM = MX + MY
      IFSX = 0
      IFSY = 0
      IFBX = 0
      IFBY = 0
      P = -ONE
C     Main loop for the different sets of knots. MPM=MX+MY is a
C     safe upper bound for the number of trials.
      DO 260 ITER = 1, MPM
         IF (NX.EQ.NMINX .AND. NY.EQ.NMINY) IER = -2
C        Find NRINTX (NRINTY) which is the number of knot intervals
C        in the X-direction (Y-direction).
         NRINTX = NX - NMINX + 1
         NRINTY = NY - NMINY + 1
C        Find NCOF, the number of B-spline coefficients for the
C        current set of knots.
         NK1X = NX - KX1
         NK1Y = NY - KY1
         NCOF = NK1X*NK1Y
C        Find the position of the additional knots which are needed
C        for the B-spline representation of S(X,Y).
         DO 160 J = 1, KX1
            TX(J) = XB
            TX(NX-J+1) = XE
  160    CONTINUE
         DO 180 J = 1, KY1
            TY(J) = YB
            TY(NY-J+1) = YE
  180    CONTINUE
C        Find the least-squares spline SINF(X,Y) and calculate for
C        each knot interval TX(J+KX)<=X<=TX(J+KX+1)
C        (TY(J+KY)<=Y<=TY(J+KY+1)) the sum of squared residuals
C        FPINTX(J),J=1,2,...,NX-2*KX-1 (FPINTY(J),J=1,2,...,NY-2*KY-1)
C        for the data points having their abscissa (ordinate)-
C        value belonging to that interval.
C        FP gives the total sum of squared residuals.
         CALL E02DCY(IFSX,IFSY,IFBX,IFBY,X,MX,Y,MY,Z,KX,KY,TX,NX,TY,NY,
     *               P,C,NC,FP,FPINTX,FPINTY,MM,MYNX,KX1,KX2,KY1,KY2,
     *               WRK(LSX),WRK(LSY),WRK(LRI),WRK(LQ),WRK(LAX),
     *               WRK(LAY),WRK(LBX),WRK(LBY),NRX,NRY)
         IF (IER.EQ.-2) FP0 = FP
C        Test whether the least-squares spline is an acceptable
C        solution.
         FPMS = FP - S
         IF (ABS(FPMS).LT.ACC) THEN
            GO TO 360
C           If F(P=INF) < S, we accept the choice of knots.
         ELSE IF (FPMS.LT.ZERO) THEN
            GO TO 300
C           If NX=NMAXX and NY=NMAXY, SINF(X,Y) is an interpolating
C           spline.
         ELSE IF (NX.EQ.NMAXX .AND. NY.EQ.NMAXY) THEN
            GO TO 280
C           Increase the number of knots.
C           If NX=NXE and NY=NYE we cannot further increase the number
C           of knots because of the storage capacity limitation.
         ELSE IF (NX.EQ.NXE .AND. NY.EQ.NYE) THEN
            IER = 4
            GO TO 360
         ELSE
            IER = 0
C           Adjust the parameter reducx or reducy according to the
C           direction in which the last added knots were located.
            IF (LASTDI.LT.0) THEN
               REDUCX = FPOLD - FP
            ELSE IF (LASTDI.GT.0) THEN
               REDUCY = FPOLD - FP
            END IF
C           Store the sum of squared residuals for the current set
C           of knots.
            FPOLD = FP
C           Find NPLX, the number of knots we should add in the
C           X-direction.
            NPLX = 1
            IF (NX.NE.NMINX) THEN
               NPL1 = NPLUSX*2
               IF (REDUCX.GT.ACC) NPL1 = NPLUSX*FPMS/REDUCX
               NPLX = MIN(NPLUSX*2,MAX(NPL1,NPLUSX/2,1))
            END IF
C           Find NPLY, the number of knots we should add in the
C           Y-direction.
            NPLY = 1
            IF (NY.NE.NMINY) THEN
               NPL1 = NPLUSY*2
               IF (REDUCY.GT.ACC) NPL1 = NPLUSY*FPMS/REDUCY
               NPLY = MIN(NPLUSY*2,MAX(NPL1,NPLUSY/2,1))
            END IF
            IF (NPLX.EQ.NPLY) THEN
               IF ((LASTDI.LT.0 .OR. NX.EQ.NXE) .AND. NY.NE.NYE)
     *             GO TO 220
            ELSE IF (NPLX.GE.NPLY .AND. NY.NE.NYE) THEN
               GO TO 220
            ELSE IF (NX.EQ.NXE .AND. NY.NE.NYE) THEN
               GO TO 220
            END IF
C           Addition in the X-direction.
            LASTDI = -1
            NPLUSX = NPLX
            IFSX = 0
            DO 200 L = 1, NPLUSX
C              Add a new knot in the X-direction
               CALL E02BEY(X,MX,TX,NX,FPINTX,NRDATX,NRINTX,NXEST)
C              Test whether we cannot further increase the number
C              of knots in the X-direction.
               IF (NX.EQ.NXE) GO TO 260
  200       CONTINUE
            GO TO 260
C           Addition in the Y-direction.
  220       LASTDI = 1
            NPLUSY = NPLY
            IFSY = 0
            DO 240 L = 1, NPLUSY
C              Add a new knot in the Y-direction.
               CALL E02BEY(Y,MY,TY,NY,FPINTY,NRDATY,NRINTY,NYEST)
C              Test whether we cannot further increase the number of
C              knots in the Y-direction.
               IF (NY.EQ.NYE) GO TO 260
  240       CONTINUE
         END IF
  260 CONTINUE
      GO TO 300
  280 IER = -1
      GO TO 360
C     Restart the computations with the new set of knots.
C     Test whether the least-squares polynomial is a solution of our
C     approximation problem.
  300 IF (IER.NE.-2) THEN
C        Part 2: Determination of the smoothing spline SP(X,Y)
C        *****************************************************
C        We have determined the number of knots and their position.
C        We now compute the B-spline coefficients of the smoothing
C        spline SP(X,Y). This smoothing spline varies with the
C        parameter P in such a way that F(P) =
C        SUMI=1,MX(SUMJ=1,MY((Z(I,J)-SP(X(I),Y(J)))**2) is a
C        continuous, strictly decreasing function of P. Moreover the
C        least-squares polynomial corresponds to P=0 and the
C        least-squares spline to P=infinity. Iteratively we then have
C        to determine the positive value of P such that F(P)=S. The
C        process which is proposed here makes use of rational
C        interpolation. F(P) is approximated by a rational function
C        R(P)=(U*P+V)/(P+W); three values of P (P1,P2,P3) with
C        corresponding values of F(P) (F1=F(P1)-S,F2=F(P2)-S,
C        F3=F(P3)-S) are used to calculate the new value of P such
C        that R(P)=S. Convergence is guaranteed by taking
C        F1 > 0 and F3 < 0.
C
C        Initial value for P.
         P1 = ZERO
         F1 = FP0 - S
         P3 = -ONE
         F3 = FPMS
         P = ONE
         ICH1 = 0
         ICH3 = 0
C        Iteration process to find the root of F(P)=S.
         DO 320 ITER = 1, MAXIT
C           Find the smoothing spline SP(X,Y) and the corresponding
C           sum of squared residuals FP.
            CALL E02DCY(IFSX,IFSY,IFBX,IFBY,X,MX,Y,MY,Z,KX,KY,TX,NX,TY,
     *                  NY,P,C,NC,FP,FPINTX,FPINTY,MM,MYNX,KX1,KX2,KY1,
     *                  KY2,WRK(LSX),WRK(LSY),WRK(LRI),WRK(LQ),WRK(LAX),
     *                  WRK(LAY),WRK(LBX),WRK(LBY),NRX,NRY)
C           Test whether the approximation SP(X,Y) is an acceptable
C           solution.
            FPMS = FP - S
            IF (ABS(FPMS).LT.ACC) THEN
               GO TO 360
C              Test whether the maximum allowable number of iterations
C              has been reached.
            ELSE IF (ITER.EQ.MAXIT) THEN
               GO TO 340
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
                     GO TO 320
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
                     GO TO 320
C                    Test whether the iteration process proceeds as
C                    theoretically expected.
                  ELSE IF (F2.GT.ZERO) THEN
                     ICH1 = 1
                  END IF
               END IF
               IF (F2.GE.F1 .OR. F2.LE.F3) THEN
                  GO TO 340
               ELSE
C                 Find the new value of P.
                  P = E02BEX(P1,F1,P2,F2,P3,F3)
               END IF
            END IF
  320    CONTINUE
  340    IER = 5
      END IF
  360 IF (IER.EQ.-1) FP = ZERO
C     The spline is an interpolant if IER = -1, so set FP = ZERO.
      IER = MAX(IER,0)
      RETURN
      END
