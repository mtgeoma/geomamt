      DOUBLE PRECISION FUNCTION E01BHY(N,X,F,D,INCFD,SKIP,IA,IB,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16A REVISED. IER-982 (JUN 1993).
C
C     ------------------------------------------------------------------
C
C     Derived from PCHIP routine
C          PCHID:  Piecewise Cubic Hermite Integrator, Data Limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
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
C     N -- (input) number of data points.  (Error return if N.lt.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .lt. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.lt.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in E01BEF).
C           SKIP will be set to .TRUE. on return with IERR = 0 or 4.
C
C     IA,IB -- (input) indices in X-array for the limits of integration.
C           Both must be in the range [1,N].  (Error return if not.)
C           No restrictions on their relative values.
C
C     IERR -- (output) error flag.
C              IERR = 0  Normal return.
C              IERR = 1  if N.lt.2 .
C              IERR = 2  if INCFD.lt.1 .
C              IERR = 3  if the X-array is not strictly increasing.
C              IERR = 4  if IA or IB is out of range.
C                (Value has not been computed if IERR .ne. 0.)
C               Note:  the above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C     ------------------------------------------------------------------
C
C     Programmed by:  Fred N. Fritsch,  FTS 532-4275,
C                     Mathematics and Statistics Division,
C                     Lawrence Livermore National Laboratory.
C
C     ------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, HALF, SIX
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0,SIX=6.0D0)
C     .. Scalar Arguments ..
      INTEGER                          IA, IB, IERR, INCFD, N
      LOGICAL                          SKIP
C     .. Array Arguments ..
      DOUBLE PRECISION                 D(INCFD,N), F(INCFD,N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 H, SUM, VALUE
      INTEGER                          I, IUP, LOW
C     .. Intrinsic Functions ..
      INTRINSIC                        MAX, MIN
C     .. Executable Statements ..
      E01BHY = 0.0D0
C     Validity-check arguments.
      IF ( .NOT. SKIP) THEN
C
         IF (N.LT.2) THEN
C
C           Error returns.
C           N.lt.2 return.
            IERR = 1
            RETURN
         ELSE IF (INCFD.LT.1) THEN
C
C           INCFD.lt.1 return.
            IERR = 2
            RETURN
         ELSE
            DO 20 I = 2, N
               IF (X(I).LE.X(I-1)) GO TO 40
   20       CONTINUE
            GO TO 60
C
   40       CONTINUE
C           X-array not strictly increasing.
            IERR = 3
            RETURN
         END IF
      END IF
C     Function definition is ok, go on.
   60 CONTINUE
      SKIP = .TRUE.
      IF ((IA.GE.1) .AND. (IA.LE.N)) THEN
         IF ((IB.GE.1) .AND. (IB.LE.N)) THEN
            IERR = 0
C           Compute integral value.
            IF (IA.EQ.IB) THEN
               VALUE = ZERO
            ELSE
               LOW = MIN(IA,IB)
               IUP = MAX(IA,IB) - 1
               SUM = ZERO
               DO 80 I = LOW, IUP
                  H = X(I+1) - X(I)
                  SUM = SUM + H*((F(1,I)+F(1,I+1))+(D(1,I)-D(1,I+1))
     *                  *(H/SIX))
   80          CONTINUE
               VALUE = HALF*SUM
               IF (IA.GT.IB) VALUE = -VALUE
            END IF
C
C           Normal return.
            E01BHY = VALUE
            RETURN
         END IF
      END IF
C
C     IA or IB out of range return.
      IERR = 4
      RETURN
      END
