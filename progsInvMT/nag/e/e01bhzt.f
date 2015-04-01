      DOUBLE PRECISION FUNCTION E01BHZ(N,X,F,D,INCFD,A,B,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ------------------------------------------------------------------
C
C     Derived from PCHIP routine
C          PCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary Limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [A, B].
C
C     To provide compatibility with E01BEZ, includes an
C     increment between successive values of the F- and D-arrays.
C
C     ------------------------------------------------------------------
C
C     Parameters:
C
C     VALUE -- (output) value of the requested integral.
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
C     A,B -- (input) the limits of integration.
C           Note:  there is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     ------------------------------------------------------------------
C
C     Programmed by:  Fred N. Fritsch,  FTS 532-4275,
C                     Mathematics and Statistics Division,
C                     Lawrence Livermore National Laboratory.
C
C     ------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO
      PARAMETER                        (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
      INTEGER                          IERR, INCFD, N
C     .. Array Arguments ..
      DOUBLE PRECISION                 D(INCFD,N), F(INCFD,N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 VALUE, XA, XB
      INTEGER                          I, IA, IB, IERD, IERV, IL, IR
      LOGICAL                          SKIP
C     .. External Functions ..
      DOUBLE PRECISION                 E01BHX, E01BHY
      EXTERNAL                         E01BHX, E01BHY
C     .. Intrinsic Functions ..
      INTRINSIC                        MAX, MIN
C     .. Executable Statements ..
C     Assume that arguments have been checked by E01BHF.
C
      IERR = 0
      IF ((A.LT.X(1)) .OR. (A.GT.X(N))) IERR = IERR - 1
      IF ((B.LT.X(1)) .OR. (B.GT.X(N))) IERR = IERR - 2
C     SKIP changed to local variable for NAG version.
      SKIP = .FALSE.
C     Compute integral value.
      IF (A.EQ.B) THEN
         VALUE = ZERO
      ELSE
         XA = MIN(A,B)
         XB = MAX(A,B)
         IF (XB.LE.X(2)) THEN
C           Interval is to left of X(2), so use first cubic.
            VALUE = E01BHX(X(1),X(2),F(1,1),F(1,2),D(1,1),D(1,2),A,B,
     *              IERV)
            IF (IERV.GT.0) GO TO 60
         ELSE IF (XA.GE.X(N-1)) THEN
C           Interval is to right of X(N-1), so use last cubic.
            VALUE = E01BHX(X(N-1),X(N),F(1,N-1),F(1,N),D(1,N-1),D(1,N),
     *              A,B,IERV)
            IF (IERV.GT.0) GO TO 60
         ELSE
C           'Normal' case -- XA.lt.XB, XA.lt.X(N-1), XB.gt.X(2).
C           Locate IA and IB such that
C           X(IA-1).lt.XA.le.X(IA).le.X(IB).le.XB.le.X(IB+1)
            IA = 1
            DO 20 I = 1, N - 1
               IF (XA.GT.X(I)) IA = I + 1
   20       CONTINUE
C           IA = 1 implies XA.lt.X(1) .  Otherwise,
C           IA is largest index such that X(IA-1).lt.XA,.
            IB = N
            DO 40 I = N, IA, -1
               IF (XB.LT.X(I)) IB = I - 1
   40       CONTINUE
C           IB = N implies XB.gt.X(N) .  Otherwise,
C           IB is smallest index such that XB.lt.X(IB+1) .
C           Compute the integral.
            IERV = 0
            IF (IB.LT.IA) THEN
C              This means IB = IA-1 and
C              [A,B] is a subset of [X(IB),X(IA)].
               VALUE = E01BHX(X(IB),X(IA),F(1,IB),F(1,IA),D(1,IB),D(1,
     *                 IA),A,B,IERV)
               IF (IERV.GT.0) GO TO 60
            ELSE
C              First compute integral over [X(IA),X(IB)].
               IF (IB.EQ.IA) THEN
                  VALUE = ZERO
               ELSE
                  VALUE = E01BHY(N,X,F,D,INCFD,SKIP,IA,IB,IERD)
                  IF (IERD.GT.0) GO TO 80
               END IF
C              Then add on integral over [XA,X(IA)].
               IF (XA.LT.X(IA)) THEN
                  IL = MAX(1,IA-1)
                  IR = IL + 1
                  VALUE = VALUE + E01BHX(X(IL),X(IR),F(1,IL),F(1,IR),
     *                    D(1,IL),D(1,IR),XA,X(IA),IERV)
                  IF (IERV.GT.0) GO TO 60
               END IF
C              Then add on integral over [X(IB),XB].
               IF (XB.GT.X(IB)) THEN
                  IR = MIN(IB+1,N)
                  IL = IR - 1
                  VALUE = VALUE + E01BHX(X(IL),X(IR),F(1,IL),F(1,IR),
     *                    D(1,IL),D(1,IR),X(IB),XB,IERV)
                  IF (IERV.GT.0) GO TO 60
               END IF
C              Finally, adjust sign if necessary.
               IF (A.GT.B) VALUE = -VALUE
            END IF
         END IF
      END IF
C
C     Normal return.
      E01BHZ = VALUE
      GO TO 100
   60 CONTINUE
C
C     Trouble in E01BHX.  (Should never occur.)
      IERR = 5
      GO TO 100
   80 CONTINUE
C
C     Trouble in E01BHY.  (Should never occur.)
      IERR = 6
  100 CONTINUE
C
      RETURN
C
      END
