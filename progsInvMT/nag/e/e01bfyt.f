      SUBROUTINE E01BFY(X1,X2,F1,F2,D1,D2,NE,XE,FE,NEXT,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ------------------------------------------------------------------
C
C     Derived from PCHIP routine
C          CHFEV:  Cubic Hermite Function Evaluator
C
C     Evaluates the cubic polynomial determined by function values
C     F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points
C     XE(J), J=1(1)NE.
C
C     ------------------------------------------------------------------
C
C     Parameters:
C
C     X1,X2 -- (input) endpoints of interval of definition of cubic.
C           (Error return if  X1.eq.X2 .)
C
C     F1,F2 -- (input) values of function at X1 and X2, respectively.
C
C     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.lt.1 .)
C
C     XE -- (input) real array of points at which the functions are to
C           be evaluated.  If any of the XE are outside the interval
C           [X1,X2], a warning error is returned.
C
C     FE -- (output) real array of values of the cubic function defined
C           by  X1,X2, F1,F2, D1,D2  at the points  XE.
C
C     NEXT -- (output) integer array indicating number of extrapolation
C           points:
C            NEXT(1) = number of evaluation points to left of interval.
C            NEXT(2) = number of evaluation points to right of interval.
C
C     IERR -- (output) error flag.
C              IERR = 0  Normal return.
C              IERR = 1  if NE.lt.1 .
C              IERR = 2  if X1.eq.X2 .
C              (Output arrays have not been changed if IERR .gt. 0)
C
C     ------------------------------------------------------------------
C
C     Programmed by:  Fred N. Fritsch,  FTS 532-4275,
C                     Mathematics and Statistics Division,
C                     Lawrence Livermore National Laboratory.
C
C     ------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D1, D2, F1, F2, X1, X2
      INTEGER           IERR, NE
C     .. Array Arguments ..
      DOUBLE PRECISION  FE(NE), XE(NE)
      INTEGER           NEXT(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2, C3, DEL1, DEL2, DELTA, H, X, XMA, XMI
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C     Validity-check arguments.
      IF (NE.LT.1) GO TO 40
      H = X2 - X1
      IF (H.EQ.ZERO) GO TO 60
C     Initialize.
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = MIN(ZERO,H)
      XMA = MAX(ZERO,H)
C     Compute cubic coefficients (expanded about X1).
      DELTA = (F2-F1)/H
      DEL1 = (D1-DELTA)/H
      DEL2 = (D2-DELTA)/H
C     (DELTA is no longer needed.)
      C2 = -(DEL1+DEL1+DEL2)
      C3 = (DEL1+DEL2)/H
C     (H, DEL1 and DEL2 are no longer needed.)
C     Evaluation loop.
      DO 20 I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1+X*(C2+X*C3))
C        Count extrapolation points.
         IF (X.LT.XMI) THEN
            NEXT(1) = NEXT(1) + 1
         ELSE IF (X.GT.XMA) THEN
            NEXT(2) = NEXT(2) + 1
         END IF
   20 CONTINUE
C
C     Normal return.
      RETURN
   40 CONTINUE
C
C     Error returns.
C     NE.lt.1 return.
      IERR = 1
      RETURN
   60 CONTINUE
C
C     X1.eq.X2 return.
      IERR = 2
      RETURN
      END
