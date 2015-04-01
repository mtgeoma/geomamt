      SUBROUTINE E04MBV(LP,N,NCTOTL,NROWA,BIGBND,FEAMIN,NUMINF,SUMINF,
     *                  ISTATE,A,BL,BU,CVEC,FEATOL,GRAD,X)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C *********************************************************************
C     IF NUMINF .GT. 0,  E04MBV  FINDS THE NUMBER AND WEIGHTED SUM OF
C     INFEASIBILITIES FOR THE BOUNDS AND LINEAR CONSTRAINTS. AN
C     APPROPRIATE GRADIENT VECTOR IS RETURNED IN  GRAD.
C     IF NUMINF = 0,  AND IF AN LP PROBLEM IS BEING SOLVED,  GRAD  WILL
C     BE LOADED WITH THE TRUE LINEAR OBJECTIVE.
C
C     POSITIVE VALUES OF  ISTATE(J)  WILL NOT BE ALTERED.  THESE MEAN
C     THE FOLLOWING...
C
C            1          2         3
C        A*X = BL   A*X = BU   BL = BU
C
C     OTHER VALUES OF  ISTATE(J)  WILL BE RESET AS FOLLOWS...
C        A*X LT BL   A*X GT BU   A*X FREE
C           - 2         - 1         0
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF SEPTEMBER 1981.  REV. OCT. 1982. JAN. 1983.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, FEAMIN, SUMINF
      INTEGER           N, NCTOTL, NROWA, NUMINF
      LOGICAL           LP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), BL(NCTOTL), BU(NCTOTL), CVEC(N),
     *                  FEATOL(NCTOTL), GRAD(N), X(N)
      INTEGER           ISTATE(NCTOTL)
C     .. Local Scalars ..
      DOUBLE PRECISION  ATX, FEASJ, S, WEIGHT, ZERO
      INTEGER           J, K, LROWA
      LOGICAL           NOLOW, NOUPP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FBF, DCOPY, DAXPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
C
      LROWA = NROWA*(N-1) + 1
      IF (NUMINF.EQ.0) GO TO 80
      NUMINF = 0
      SUMINF = ZERO
      CALL F06FBF(N,ZERO,GRAD,1)
C
      DO 60 J = 1, NCTOTL
C
C        DO NOTHING IF THE VARIABLE OR CONSTRAINT IS AT A BOUND.
C
         IF (ISTATE(J).GT.0) GO TO 60
         FEASJ = FEATOL(J)
         NOLOW = BL(J) .LE. (-BIGBND)
         NOUPP = BU(J) .GE. BIGBND
         K = J - N
         IF (J.LE.N) ATX = X(J)
         IF (J.GT.N) ATX = DDOT(N,A(K,1),NROWA,X,1)
         ISTATE(J) = 0
C
C        SEE IF THE LOWER BOUND IS VIOLATED.
C
         IF (NOLOW) GO TO 20
         S = BL(J) - ATX
         IF (S.LE.FEASJ) GO TO 20
         ISTATE(J) = -2
         WEIGHT = -FEAMIN/FEASJ
         GO TO 40
C
C        SEE IF THE UPPER BOUND IS VIOLATED.
C
   20    IF (NOUPP) GO TO 60
         S = ATX - BU(J)
         IF (S.LE.FEASJ) GO TO 60
         ISTATE(J) = -1
         WEIGHT = FEAMIN/FEASJ
C
C        ADD THE INFEASIBILITY.
C
   40    NUMINF = NUMINF + 1
         SUMINF = SUMINF + ABS(WEIGHT)*S
         IF (J.LE.N) GRAD(J) = WEIGHT
         IF (J.GT.N) CALL DAXPY(N,WEIGHT,A(K,1),NROWA,GRAD,1)
   60 CONTINUE
C
C     IF FEASIBLE, INSTALL TRUE OBJECTIVE.
C
   80 IF (LP .AND. NUMINF.EQ.0) CALL DCOPY(N,CVEC,1,GRAD,1)
      RETURN
C
C     END OF E04MBV  ( LPGRAD )
      END
