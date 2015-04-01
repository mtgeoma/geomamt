      SUBROUTINE E01BGZ(N,X,F,D,INCFD,M,PX,PF,PD,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ------------------------------------------------------------------
C
C     Derived from PCHIP routine
C          PCHFD:  Piecewise Cubic Hermite Function and Derivative
C                  Evaluator
C
C     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
C     gether with its first derivative, at the points  PX(J), J=1(1)M.
C
C     If only function values are required, use E01BFZ instead.
C
C     To provide compatibility with E01BEZ, includes an
C     increment between successive values of the F- and D-arrays.
C
C     ------------------------------------------------------------------
C
C     Parameters:
C
C     N -- (input) number of data points.  N.ge.2.
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .lt. X(I),  I = 2(1)N.
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           INCFD.ge.1.
C
C     M -- (input) number of evaluation points.  M.ge.1.
C
C     PX -- (input) real array of points at which the functions are to
C           be evaluated.
C           Notes:
C           1. The evaluation will be most efficient if the elements
C              of PX are increasing relative to X;
C              that is,   PX(J) .ge. X(I)
C              implies    PX(K) .ge. X(I),  all K.GE.J .
C           2. If any of the PX are outside the interval [X(1),X(N)],
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C     PF -- (output) real array of values of the cubic hermite function
C           defined by  N, X, F, D  at the points  PX.
C
C     PD -- (output) real array of values of the first derivative of
C           the same function at the points  PX.
C
C     ------------------------------------------------------------------
C
C     Programmed by:  Fred N. Fritsch,  FTS 532-4275,
C                     Mathematics and Statistics Division,
C                     Lawrence Livermore National Laboratory.
C
C     ------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IERR, INCFD, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(INCFD,N), F(INCFD,N), PD(M), PF(M), PX(M),
     *                  X(N)
C     .. Local Scalars ..
      INTEGER           I, IERC, IR, J, JFIRST, NJ
C     .. Local Arrays ..
      INTEGER           NEXT(2)
C     .. External Subroutines ..
      EXTERNAL          E01BGY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C     Assume that arguments have been checked by E01BGF.
C
C     Loop over intervals.  (   Interval index is  IL = IR-1  . )
C                           ( Interval is X(IL).le.X.lt.X(IR) . )
      JFIRST = 1
      IR = 2
   20 CONTINUE
C     Skip out of loop if have processed all evaluation points.
      IF (JFIRST.GT.M) THEN
         GO TO 240
      ELSE
C        Locate all points in interval.
         DO 40 J = JFIRST, M
            IF (PX(J).GE.X(IR)) GO TO 60
   40    CONTINUE
         J = M + 1
         GO TO 80
C        Have located first point beyond interval.
   60    CONTINUE
         IF (IR.EQ.N) J = M + 1
C
   80    CONTINUE
         NJ = J - JFIRST
C        Skip evaluation if no points in interval.
         IF (NJ.NE.0) THEN
C           Evaluate cubic at PX(I),  I = JFIRST (1) J-1 .
            CALL E01BGY(X(IR-1),X(IR),F(1,IR-1),F(1,IR),D(1,IR-1),D(1,
     *                  IR),NJ,PX(JFIRST),PF(JFIRST),PD(JFIRST),NEXT,
     *                  IERC)
            IF (IERC.GT.0) THEN
               GO TO 220
            ELSE
C
               IF (NEXT(2).GT.0) THEN
C                 In the current set of PX-points, There Are NEXT(2)
C                 To The Right Of X(IR).
                  IF (IR.EQ.N) THEN
C                    These are actually extrapolation points.
                     IERR = IERR - NEXT(2)
                  ELSE
                     GO TO 180
                  END IF
               END IF
C
               IF (NEXT(1).GT.0) THEN
C                 In the current set of PX-points, there are NEXT(1)
C                 to the left of X(IR-1).
                  IF (IR.EQ.2) THEN
C                    These are actually extrapolation points.
                     IERR = IERR - NEXT(1)
                  ELSE
C                    PX is not ordered relative to X, so must adjust
C                    evaluation interval.
C                    First, locate first point to left of X(IR-1).
                     DO 100 I = JFIRST, J - 1
                        IF (PX(I).LT.X(IR-1)) GO TO 120
  100                CONTINUE
                     GO TO 200
C
  120                CONTINUE
C                    Reset J.  (This will be the new JFIRST.)
                     J = I
C                    Now find out how far to back up in the X-array.
                     DO 140 I = 1, IR - 1
                        IF (PX(J).LT.X(I)) GO TO 160
  140                CONTINUE
C                    NB:  can never drop through here,
C                         since PX(J).lt.X(IR-1).
  160                CONTINUE
C                    At this point, either  PX(J) .lt. X(1)
C                    or      X(I-1) .le. PX(J) .lt. X(I) .
C                    Reset IR, recognizing that it will be
C                    incremented before cycling.
                     IR = MAX(1,I-1)
                  END IF
               END IF
C
               JFIRST = J
            END IF
         END IF
C        End of IR-loop.
         IR = IR + 1
         IF (IR.LE.N) THEN
            GO TO 20
         ELSE
            GO TO 240
         END IF
      END IF
C     We should never have gotten here.
  180 CONTINUE
      GO TO 220
C     Note:  cannot drop through here unless there is an error
C            in E01BGY.
  200 CONTINUE
C
  220 CONTINUE
C     Error return from E01BGY.
C     *** This case should never occur ***
      IERR = 6
C
C     Normal return.
  240 CONTINUE
      RETURN
C
      END
