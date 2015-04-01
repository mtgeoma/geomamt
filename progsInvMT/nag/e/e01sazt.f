      SUBROUTINE E01SAZ(N,X,Y,Z,IADJ,IEND,EPS,NIT,ZXZY,COINC1,COINC2,
     *                  IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: GRADG.
C     Given a triangulation of N nodes in the plane with
C     associated data values, this routine uses a global method
C     to compute estimated gradients at the nodes.  The method
C     consists of minimizing a quadratic functional Q(G) over
C     the N-vector G of gradients where Q approximates the lin-
C     earized curvature of an interpolant F over the triangula-
C     tion.  The restriction of F to an arc of the triangulation
C     is taken to be the Hermite cubic interpolant of the data
C     values and tangential gradient components at the end-
C     points of the arc, and Q is the sum of the linearized
C     curvatures of F along the arcs -- the integrals over the
C     arcs of d2F(T)**2 where d2F(t) is the second derivative
C     of F with respect to distance T along the arc.  This min-
C     imization problem corresponds to an order 2N symmetric
C     positive-definite sparse linear system which is solved for
C     the X and Y partial derivatives by the block Gauss-Seidel
C     method with 2 by 2 blocks.
C
C     Input Parameters - N - number of nodes.  N .ge. 3.
C
C                  X,Y - cartesian coordinates of the nodes.
C
C                    Z - data values at the nodes.  Z(I) is
C                        associated with (X(I),Y(I)).
C
C            IADJ,IEND - data structure defining the trian-
C                        gulation.  See subroutine E01SAY.
C
C                  EPS - nonnegative convergence criterion.
C                        The method is terminated when the
C                        maximum change in a gradient compo-
C                        nent between iterations is at most
C                        EPS.  EPS = 1.E-2 is sufficient for
C                        effective convergence.
C
C     The above parameters are not altered by this routine.
C
C                  NIT - maximum number of Gauss-Seidel
C                        iterations to be applied.  This
C                        maximum will likely be achieved if
C                        EPS is smaller than the machine
C                        precision.  Optimal efficiency was
C                        achieved in testing with EPS = 0
C                        and NIT = 3 or 4.
C
C                 ZXZY - 2 by N array whose columns contain
C                        initial estimates of the partial
C                        derivatives (zero vectors are
C                        sufficient).
C
C     Output Parameters - NIT - number of Gauss-Seidel itera-
C                           tions employed.
C
C                    ZXZY - estimated X and Y partial deriv-
C                           atives at the nodes with X par-
C                           tials in the first row.  ZXZY is
C                           not changed if IER = 2.
C
C          COINC1, COINC2 - see description for IER = 3.
C
C                     IER - error indicator
C                           IER = 0 if the convergence cri-
C                                   terion was achieved.
C                           IER = 1 if convergence was not
C                                   achieved within nit
C                                   iterations.
C                           IER = 2 if N or EPS is out of
C                                   range or NIT .lt. 0 on
C                                   input.
C                           IER = 3 if two nodes are coincidental.
C                                   COINC1 and COINC2 are the
C                                   indices of the nodes.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NN =          local copy of N
C     MAXIT =       input value of NIT
C     ITER =        number of iterations used
C     K =           do-loop and node index
C     INDF,INDL =   IADJ indices of the first and last neighbors
C                 of K
C     INDX =        IADJ index in the range INDF,...,INDL
C     NB =          neighbor of K
C     TOL =         local copy of EPS
C     DGMAX =       maximum change in a gradient component be-
C                 tween iterations
C     NORMG =       absolute value of largest ZXZY element.
C     XK,YK,ZK =    X(K), Y(K), Z(K)
C     ZXK,ZYK =     initial values of ZXZY(1,K) and ZXZY(2,K)
C     A11,A12,A22 = matrix components of the 2 by 2 block A*DG
C                 = R where A is symmetric, DG = (DZX,DZY)
C                 is the change in the gradient at K, and R
C                 is the residual
C     R1,R2 =       components of the residual -- derivatives of
C                 Q with respect to the components of the
C                 gradient at node K
C     DELX,DELY =   components of the arc NB-K
C     DELXS,DELYS = DELX**2, DELY**2
C     D =           the distance between K and NB
C     DCUB =        D**3
C     T =           factor of R1 and R2
C     DZX,DZY =     solution of the 2 by 2 system -- change in
C                 derivatives at K from the previous iterate
C
C     .. Parameters ..
      DOUBLE PRECISION  ONEHLF, TWO, ZERO
      PARAMETER         (ONEHLF=1.5D0,TWO=2.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           COINC1, COINC2, IER, N, NIT
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N), Z(N), ZXZY(2,N)
      INTEGER           IADJ(6*N), IEND(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A11, A12, A22, D, DCUB, DELX, DELXS, DELY,
     *                  DELYS, DGMAX, DZX, DZY, NORMG, R1, R2, T, TOL,
     *                  XK, YK, ZK, ZXK, ZYK
      INTEGER           INDF, INDL, INDX, ITER, K, MAXIT, NB, NN
C     .. External Functions ..
      DOUBLE PRECISION  F06BNF
      EXTERNAL          F06BNF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      NN = N
      TOL = EPS
      MAXIT = NIT
      COINC1 = 0
      COINC2 = 0
      IF (NN.LT.3 .OR. TOL.LT.ZERO .OR. MAXIT.LT.0) THEN
C        Parameter out of range
         NIT = 0
         IER = 2
      ELSE
         NORMG = ZERO
         DO 20 K = 1, NN
            NORMG = MAX(NORMG,ABS(ZXZY(1,K)),ABS(ZXZY(2,K)))
   20    CONTINUE
         ITER = 0
C        Top of iteration loop
   40    IF (ITER.GT.MAXIT) THEN
C           Method has failed to converge within MAXIT iterations.
            IER = 0
         ELSE
            DGMAX = ZERO
            INDL = 0
            DO 80 K = 1, NN
               XK = X(K)
               YK = Y(K)
               ZK = Z(K)
               ZXK = ZXZY(1,K)
               ZYK = ZXZY(2,K)
C              Initialize components of the 2 by 2 system
               A11 = ZERO
               A12 = ZERO
               A22 = ZERO
               R1 = ZERO
               R2 = ZERO
C              Loop on neighbors NB of K
               INDF = INDL + 1
               INDL = IEND(K)
               DO 60 INDX = INDF, INDL
                  NB = IADJ(INDX)
                  IF (NB.NE.0) THEN
C                    Compute the components of arc NB-K
                     DELX = X(NB) - XK
                     DELY = Y(NB) - YK
                     DELXS = DELX*DELX
                     DELYS = DELY*DELY
                     D = F06BNF(DELX,DELY)
                     IF (D.EQ.ZERO) THEN
C                       Nodes K and NB are coincidental.
                        IER = 3
                        COINC1 = K
                        COINC2 = NB
                        RETURN
                     ELSE
                        DCUB = D*D*D
C                       Update the system components for node NB.
                        A11 = A11 + DELXS/DCUB
                        A12 = A12 + DELX*DELY/DCUB
                        A22 = A22 + DELYS/DCUB
                        T = (ONEHLF*(Z(NB)-ZK)-((ZXZY(1,NB)/TWO+ZXK)
     *                      *DELX+(ZXZY(2,NB)/TWO+ZYK)*DELY))/DCUB
                        R1 = R1 + T*DELX
                        R2 = R2 + T*DELY
                     END IF
                  END IF
   60          CONTINUE
C              Solve the 2 by 2 system and update DGMAX
               DZY = (A11*R2-A12*R1)/(A11*A22-A12*A12)
               DZX = (R1-A12*DZY)/A11
               DGMAX = MAX(DGMAX,ABS(DZX),ABS(DZY))
C              Update the partials at node K
               ZXZY(1,K) = ZXK + DZX
               ZXZY(2,K) = ZYK + DZY
   80       CONTINUE
C           Increment ITER and test for convergence
            ITER = ITER + 1
            IF (DGMAX.GT.TOL+TOL*NORMG) GO TO 40
C           Method has converged
            IER = 0
         END IF
      END IF
      NIT = ITER
      RETURN
      END
