      SUBROUTINE F08JEX(A,B,C,RT1,RT2)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAE2(A,B,C,RT1,RT2)
C
C  Purpose
C  =======
C
C  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
C     [  A   B  ]
C     [  B   C  ].
C  On return, RT1 is the eigenvalue of larger absolute value, and RT2
C  is the eigenvalue of smaller absolute value.
C
C  Arguments
C  =========
C
C  A       (input) DOUBLE PRECISION
C          The (1,1) entry of the 2-by-2 matrix.
C
C  B       (input) DOUBLE PRECISION
C          The (1,2) and (2,1) entries of the 2-by-2 matrix.
C
C  C       (input) DOUBLE PRECISION
C          The (2,2) entry of the 2-by-2 matrix.
C
C  RT1     (output) DOUBLE PRECISION
C          The eigenvalue of larger absolute value.
C
C  RT2     (output) DOUBLE PRECISION
C          The eigenvalue of smaller absolute value.
C
C  Further Details
C  ===============
C
C  RT1 is accurate to a few ulps barring over/underflow.
C
C  RT2 may be inaccurate if there is massive cancellation in the
C  determinant A*C-B*B; higher precision or correctly rounded or
C  correctly truncated arithmetic would be needed to compute RT2
C  accurately in all cases.
C
C  Overflow is possible only if RT1 is within a factor of 5 of overflow.
C  Underflow is harmless if the input data is 0 or exceeds
C     underflow_threshold / macheps.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D0)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, C, RT1, RT2
C     .. Local Scalars ..
      DOUBLE PRECISION  AB, ACMN, ACMX, ADF, DF, RT, SM, TB
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
C     Compute the eigenvalues
C
      SM = A + C
      DF = A - C
      ADF = ABS(DF)
      TB = B + B
      AB = ABS(TB)
      IF (ABS(A).GT.ABS(C)) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF (ADF.GT.AB) THEN
         RT = ADF*SQRT(ONE+(AB/ADF)**2)
      ELSE IF (ADF.LT.AB) THEN
         RT = AB*SQRT(ONE+(ADF/AB)**2)
      ELSE
C
C        Includes case AB=ADF=0
C
         RT = AB*SQRT(TWO)
      END IF
      IF (SM.LT.ZERO) THEN
         RT1 = HALF*(SM-RT)
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE IF (SM.GT.ZERO) THEN
         RT1 = HALF*(SM+RT)
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE
C
C        Includes case RT1 = RT2 = 0
C
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
C
C     End of F08JEX (DLAE2)
C
      END
