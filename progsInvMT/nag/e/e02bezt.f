      SUBROUTINE E02BEZ(IOPT,X,Y,W,M,XB,XE,K,S,NEST,TOL,MAXIT,K1,K2,N,T,
     *                  C,FP,FPINT,Z,A,B,G,Q,NRDATA,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 17 REVISED. IER-1680 (JUN 1995).
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FP, S, TOL, XB, XE
      INTEGER           IER, IOPT, K, K1, K2, M, MAXIT, N, NEST
C     .. Array Arguments ..
      DOUBLE PRECISION  A(K1,NEST+K1), B(NEST,K2), C(NEST), FPINT(NEST),
     *                  G(K2,NEST+K2), Q(M,K1), T(NEST), W(M), X(M),
     *                  Y(M), Z(NEST)
      INTEGER           NRDATA(NEST)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, COS, F1, F2, F3, FP0, FPART, FPMS, FPOLD,
     *                  P, P1, P2, P3, PINV, PIV, SIN, STORE, TERM, WI,
     *                  XI
      INTEGER           I, I2, ICH1, ICH3, IMP, IT, ITER, J, K3, L, MK1,
     *                  N8, NEW, NK1, NMAX, NMIN, NPL1, NPLUS, NRINT
      LOGICAL           RL
C     .. Local Arrays ..
      DOUBLE PRECISION  H(7), YI(1)
C     .. External Functions ..
      DOUBLE PRECISION  E02BEX
      EXTERNAL          E02BEX
C     .. External Subroutines ..
      EXTERNAL          E02BEV, E02BEW, E02BEY, DROTG, DROT, DTBSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C     Part 1: Determination of the number of knots and their position
C     ***************************************************************
C     Given a set of knots we compute the least-squares spline SINF(X).
C     If the sum F(P=INF) <= S we accept the choice of knots,
C     otherwise we have to increase their number.
C     The initial choice of knots depends on the value of S and IOPT.
C     if S=0 we have spline interpolation; in that case the number of
C     knots equals NMAX = M+K+1.
C     if S > 0 and
C      IOPT=0 we first compute the least-squares polynomial of
C      degree K; N = NMIN = 2*K+2
C      IOPT=1 we start with the set of knots found at the last
C      call of the routine, except for the case that S > FP0; then
C      we compute directly the least-squares polynomial of degree K.
C
C     Calculation of ACC, the absolute tolerance for the root of F(P)=S.
      ACC = TOL*S
C     Determine NMIN, the number of knots for polynomial approximation.
      NMIN = 2*K1
C     Determine NMAX, the number of knots for spline interpolation.
      NMAX = M + K1
      RL = .FALSE.
      IF (S.EQ.ZERO) THEN
C        If S=0, S(X) is an interpolating spline.
C        Test whether the required storage space exceeds the available
C        one.
         IF (IOPT.NE.0) FP0 = FPINT(N)
         N = NMAX
         IF (NMAX.GT.NEST) THEN
            IER = 4
            GO TO 620
         ELSE
            FPOLD = ZERO
            NPLUS = 0
         END IF
      END IF
   20 CONTINUE
C     Find the position of the interior knots in case of interpolation.
      IF (S.EQ.ZERO .OR. RL) THEN
         MK1 = M - K1
         IF (MK1.NE.0) THEN
            K3 = K/2
            IF (K3*2.EQ.K) THEN
               DO 40 L = 1, MK1
                  T(K2+L-1) = (X(K3+L-1)+X(K3+L))*HALF
   40          CONTINUE
            ELSE
               DO 60 L = 1, MK1
                  T(K2+L-1) = X(K3+L+1)
   60          CONTINUE
            END IF
         END IF
      ELSE
C        If S>0 our initial choice of knots depends on the value of
C        IOPT. If IOPT=0 or IOPT=1 and S>=FP0, we start computing the
C        least-squares polynomial of degree K which is a spline without
C        interior knots. If IOPT=1 and FP0>S we start computing the
C        least squares spline according to the set of knots found at
C        the last call of the routine.
         IF (IOPT.NE.0) THEN
            IF (N.NE.NMIN) THEN
               FP0 = FPINT(N)
               FPOLD = FPINT(N-1)
               NPLUS = NRDATA(N)
               IF (FP0.GT.S) GO TO 80
            END IF
         END IF
         N = NMIN
         FPOLD = ZERO
         NPLUS = 0
         NRDATA(1) = M - 2
      END IF
C     Main loop for the different sets of knots. M is a safe upper bound
C     for the number of trials.
   80 DO 320 ITER = 1, M
         IF (N.EQ.NMIN) IER = -2
C        Find NRINT, tne number of knot intervals.
         NRINT = N - NMIN + 1
C        Find the position of the additional knots which are needed for
C        the B-Spline representation of S(X).
         NK1 = N - K1
         DO 100 J = 1, K1
            T(J) = XB
            T(N-J+1) = XE
  100    CONTINUE
C        Compute the B-Spline coefficients of the least-squares spline
C        SINF(X). The observation matrix A is built up row by row and
C        reduced to upper triangular form by Givens transformations.
C        At the same time FP=F(P=INF) is computed.
         FP = ZERO
C        Initialize the observation matrix A.
         DO 120 I = 1, NK1
            Z(I) = ZERO
  120    CONTINUE
         DO 160 I = 1, K1
            DO 140 J = K1 - I + 1, K1 - I + NK1
               A(I,J) = ZERO
  140       CONTINUE
  160    CONTINUE
         L = K1
         DO 220 IT = 1, M
C           Fetch the current data point X(IT),Y(IT).
            XI = X(IT)
            WI = W(IT)
            YI(1) = Y(IT)*WI
C           Search for knot interval T(L) <= XI < T(L+1).
            IF (XI.GE.T(L+1) .AND. L.NE.NK1) L = L + 1
C           Evaluate the (K+1) non-zero B-Splines at XI and store them
C           in Q.
            CALL E02BEV(T,N,K,XI,L,H)
            DO 180 I = 1, K1
               Q(IT,I) = H(I)
               H(I) = H(I)*WI
  180       CONTINUE
C           Rotate the new row of the observation matrix into triangle.
            J = L - K1
            DO 200 I = 1, K1
               J = J + 1
               PIV = H(I)
               IF (PIV.NE.ZERO) THEN
C                 Calculate the parameters of the Givens transformation.
                  CALL DROTG(A(K1,J),PIV,COS,SIN)
C                 Transformations to right hand side.
                  CALL DROT(1,YI(1),1,Z(J),1,COS,-SIN)
C                 Transformations to left hand side.
                  CALL DROT(K1-I,H(I+1),1,A(K1-1,J+1),K1-1,COS,-SIN)
               END IF
  200       CONTINUE
C           Add contribution of this row to the sum of squares of
C           residual right hand sides.
            FP = FP + YI(1)*YI(1)
  220    CONTINUE
         IF (IER.EQ.-2) FP0 = FP
         FPINT(N) = FP0
         FPINT(N-1) = FPOLD
         NRDATA(N) = NPLUS
         DO 240 IMP = 1, NK1
            C(IMP) = Z(IMP)
  240    CONTINUE
C        Backward substitution to obtain the B-Spline coefficients.
         CALL DTBSV('U','N','N',NK1,K1-1,A,K1,C,1)
C        Test whether the approximation SINF(X) is an acceptable
C        solution.
         FPMS = FP - S
         IF (ABS(FPMS).LT.ACC) THEN
            GO TO 620
C           If F(P=INF) < S accept the choice of knots.
         ELSE IF (FPMS.LT.ZERO) THEN
            GO TO 360
C           If N = NMAX, SINF(X) is an interpolating spline.
         ELSE IF (N.EQ.NMAX) THEN
            GO TO 340
C           Increase the number of knots.
C           If N=NEST we cannot increase the number of knots because of
C           the storage capacity limitation.
         ELSE IF (N.EQ.NEST) THEN
            IER = 4
            GO TO 620
         ELSE
C           Determine the number of knots NPLUS we are going to add.
            IF (IER.NE.0) THEN
               NPLUS = 1
               IER = 0
            ELSE
               NPL1 = NPLUS*2
               IF (FPOLD-FP.GT.ACC) NPL1 = NPLUS*FPMS/(FPOLD-FP)
               NPLUS = MIN(NPLUS*2,MAX(NPL1,NPLUS/2,1))
            END IF
            FPOLD = FP
C           Compute the sum((W(I)*(Y(I)-S(X(I))))**2) for each knot
C           interval T(J+K) <= X(I) <= T(J+K+1) and store it in
C           FPINT(J),J=1,2,...NRINT.
            FPART = ZERO
            I = 1
            L = K2
            NEW = 0
            DO 280 IT = 1, M
               IF (X(IT).GE.T(L) .AND. L.LE.NK1) THEN
                  NEW = 1
                  L = L + 1
               END IF
               TERM = ZERO
               DO 260 J = 1, K1
                  TERM = TERM + C(L-K2+J)*Q(IT,J)
  260          CONTINUE
               TERM = (W(IT)*(TERM-Y(IT)))**2
               FPART = FPART + TERM
               IF (NEW.NE.0) THEN
                  STORE = TERM*HALF
                  FPINT(I) = FPART - STORE
                  I = I + 1
                  FPART = STORE
                  NEW = 0
               END IF
  280       CONTINUE
            FPINT(NRINT) = FPART
            DO 300 L = 1, NPLUS
C              Add a new knot.
               CALL E02BEY(X,M,T,N,FPINT,NRDATA,NRINT,NEST)
C              If N=NMAX we locate the knots as for interpolation.
               IF (N.EQ.NMAX) THEN
                  RL = .TRUE.
                  GO TO 20
C                 Test whether we cannot further increase the
C                 number of knots.
               ELSE IF (N.EQ.NEST) THEN
                  GO TO 320
               END IF
  300       CONTINUE
         END IF
  320 CONTINUE
      GO TO 360
  340 IER = -1
      GO TO 620
C     Restart the computations with the new set of knots.
C     Test whether the least-squares Kth degree polynomial is a solution
C     of our approximation problem.
  360 IF (IER.NE.-2) THEN
C        Part 2: Determination of the smoothing spline SP(X).
C        ***************************************************
C        We have determined the number of knots and their position.
C        We now compute the B-Spline coefficients of the smoothing
C        spline SP(X). The observation matrix A is extended by the
C        rows of matrix B expressing that the Kth derivative
C        discontinuities of SP(X) at the interior knots
C        T(K+2),...T(N-K-1) must be zero. The corresponding weights
C        of these additional rows are set to 1/P. Iteratively we
C        then have to determine the value of P such that
C        F(P)=SUM((W(I)*(Y(I)-SP(X(I))))**2) BE = S. We already know
C        that the least-squares Kth degree polynomial corresponds to
C        P=0, and that the least-squares spline corresponds to
C        P=infinity. The iteration process which is proposed here,
C        makes use of rational interpolation. since F(P) is a convex
C        and strictly decreasing function of P, it can be
C        approximated by a rational function R(P) = (U*P+V)/(P+W).
C        Three values of P(P1,P2,P3) with corresponding values of
C        F(P) (F1=F(P1)-S,F2=F(P2)-S,F3=F(P3)-S) are used to
C        calculate the new value of P such that R(P)=S. Convergence
C        is guaranteed by taking F1>0 AND F3<0.
C
C        Evaluate the discontinuity jump of the Kth derivative of the
C        B-Splines at the knots T(L),L=K+2,...N-K-1 and store in B.
         CALL E02BEW(T,N,K2,B,NEST)
C        Initial value for P.
         P1 = ZERO
         F1 = FP0 - S
         P3 = -ONE
         F3 = FPMS
         P = ZERO
         DO 380 I = 1, NK1
            P = P + A(K1,I)
  380    CONTINUE
         P = NK1/P
         ICH1 = 0
         ICH3 = 0
         N8 = N - NMIN
C        Iteration process to find the root of F(P) = S.
         DO 580 ITER = 1, MAXIT
C           The rows of matrix B with weight 1/P are rotated into the
C           triangularised observation matrix A which is stored in G.
            PINV = ONE/P
            DO 400 I = 1, NK1
               C(I) = Z(I)
               G(1,I+K1) = ZERO
  400       CONTINUE
            DO 440 I = 1, K1
               DO 420 J = K1 - I + 1, K1 - I + NK1
                  G(I+1,J) = A(I,J)
  420          CONTINUE
  440       CONTINUE
            DO 520 IT = 1, N8
C              The row of matrix B is rotated into triangle by Givens
C              transformation
               DO 460 I = 1, K2
                  H(I) = B(IT,I)*PINV
  460          CONTINUE
               YI(1) = ZERO
               DO 500 J = IT, NK1
                  PIV = H(1)
C                 Calculate the parameters of the Givens transformation.
                  CALL DROTG(G(K1+1,J),PIV,COS,SIN)
C                 Transformations to right hand side.
                  CALL DROT(1,YI(1),1,C(J),1,COS,-SIN)
                  IF (J.NE.NK1) THEN
                     I2 = K1
                     IF (J.GT.N8) I2 = NK1 - J
C                    Transformations to left hand side.
                     CALL DROT(I2,H(2),1,G(K1,J+1),K2-1,COS,-SIN)
                     DO 480 I = 1, I2
                        H(I) = H(I+1)
  480                CONTINUE
                     H(I2+1) = ZERO
                  END IF
  500          CONTINUE
  520       CONTINUE
C           Backward substitution to obtain the B-Spline coefficients.
            CALL DTBSV('U','N','N',NK1,K2-1,G,K2,C,1)
C           Computation of F(P).
            FP = ZERO
            L = K2
            DO 560 IT = 1, M
               IF (X(IT).GE.T(L) .AND. L.LE.NK1) L = L + 1
               TERM = ZERO
               DO 540 J = 1, K1
                  TERM = TERM + C(L-K2+J)*Q(IT,J)
  540          CONTINUE
               FP = FP + (W(IT)*(TERM-Y(IT)))**2
  560       CONTINUE
C           Test whether the approximation SP(X) is an acceptable
C           solution.
            FPMS = FP - S
            IF (ABS(FPMS).LT.ACC) THEN
               GO TO 620
C              Test whether the maximum number of iterations is reached.
            ELSE IF (ITER.EQ.MAXIT) THEN
               GO TO 600
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
                     GO TO 580
                  ELSE IF (F2.LT.ZERO) THEN
                     ICH3 = 1
                  END IF
               END IF
               IF (ICH1.EQ.0) THEN
                  IF ((F1-F2).LE.ACC) THEN
C                    Our initial choice of p is too small
                     P1 = P2
                     F1 = F2
                     P = P*25
                     IF (P3.GE.ZERO .AND. P.GE.P3) P = (P2+P3*9)/10
                     GO TO 580
                  ELSE IF (F2.GT.ZERO) THEN
                     ICH1 = 1
                  END IF
               END IF
C              Test whether the iteration process proceeds as
C              theoretically expected.
               IF (F2.GE.F1 .OR. F2.LE.F3) THEN
                  GO TO 600
               ELSE
C                 Find the new value for P.
                  P = E02BEX(P1,F1,P2,F2,P3,F3)
               END IF
            END IF
  580    CONTINUE
  600    IER = 5
      END IF
  620 IER = MAX(IER,0)
      RETURN
      END
