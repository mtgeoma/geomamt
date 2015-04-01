      SUBROUTINE E02BEY(X,M,T,N,FPINT,NRDATA,NRINT,NEST)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     Subroutine E02BEY locates an additional knot for a spline
C     of degree K and adjusts the corresponding parameters,i.e.
C     T     : the position of the knots.
C     N     : the number of knots.
C     NRINT : the number of knot intervals.
C     FPINT : the sum of squares of residual right hand sides
C             for each knot interval.
C     NRDATA: the number of data points inside each knot interval.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           M, N, NEST, NRINT
C     .. Array Arguments ..
      DOUBLE PRECISION  FPINT(NEST), T(NEST), X(M)
      INTEGER           NRDATA(NEST)
C     .. Local Scalars ..
      DOUBLE PRECISION  FPMAX
      INTEGER           IHALF, J, JBEGIN, JPOINT, K, MAXBEG, MAXPT,
     *                  NEXT, NRX, NUMBER
C     .. Executable Statements ..
C
      K = (N-NRINT-1)/2
C     Search for knot interval T(NUMBER+K) <= X <= T(NUMBER+K+1) where
C     FPINT(NUMBER) is maximal on the condition that NRDATA(NUMBER)
C     does not equal zero.
      FPMAX = ZERO
      JBEGIN = 1
      DO 20 J = 1, NRINT
         JPOINT = NRDATA(J)
         IF (FPMAX.LT.FPINT(J) .AND. JPOINT.NE.0) THEN
            FPMAX = FPINT(J)
            NUMBER = J
            MAXPT = JPOINT
            MAXBEG = JBEGIN
         END IF
         JBEGIN = JBEGIN + JPOINT + 1
   20 CONTINUE
C     Let coincide the new knot T(NUMBER+K+1) with a data point X(NRX)
C     inside the old knot interval T(NUMBER+K) <= X <= T(NUMBER+K+1).
      IHALF = MAXPT/2 + 1
      NRX = MAXBEG + IHALF
      NEXT = NUMBER + 1
C     Adjust the different parameters.
      DO 40 J = NRINT, NEXT, -1
         FPINT(J+1) = FPINT(J)
         NRDATA(J+1) = NRDATA(J)
         T(J+K+1) = T(J+K)
   40 CONTINUE
      NRDATA(NUMBER) = IHALF - 1
      NRDATA(NEXT) = MAXPT - IHALF
      FPINT(NUMBER) = (FPMAX*NRDATA(NUMBER))/MAXPT
      FPINT(NEXT) = (FPMAX*NRDATA(NEXT))/MAXPT
      T(NEXT+K) = X(NRX)
      N = N + 1
      NRINT = NRINT + 1
      RETURN
      END
