      LOGICAL FUNCTION H02BBU(X,N,INTVAR,TOLIV)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     FUNCTION DETERMINES WHETHER A SOLUTION FROM A RELAXED LP IS
C     ALSO A FEASIBLE POINT OF THE MIXED INTEGER PROGRAM.
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        TOLIV
      INTEGER                 N
C     .. Array Arguments ..
      DOUBLE PRECISION        X(*)
      INTEGER                 INTVAR(*)
C     .. Local Scalars ..
      DOUBLE PRECISION        DIFF
      INTEGER                 I
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, NINT
C     .. Executable Statements ..
C
      H02BBU = .TRUE.
      DO 20 I = 1, N
         IF (INTVAR(I).EQ.1) THEN
C            IF (ABS(X(I)-NINT(X(I))).GT.TOLIV*(1.D0+ABS(X(I)))) THEN
            DIFF = ABS(X(I)-NINT(X(I)))
            IF (DIFF.GT.TOLIV) THEN
               H02BBU = .FALSE.
               RETURN
            END IF
         END IF
   20 CONTINUE
      RETURN
      END
