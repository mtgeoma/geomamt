      SUBROUTINE F07MDX(A,B,C,RT1,RT2,CS1,SN1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             DLAEV2(A,B,C,RT1,RT2,CS1,SN1)
C
C  Purpose
C  =======
C
C  DLAEV2 computes the eigendecomposition of a 2 by 2 symmetric matrix:
C     ((A,B);(B,C))
C
C  RT1 is the eigenvalue of larger absolute value, RT2 is the eigenvalue
C  of smaller absolute value, and (CS1,SN1) is the unit right
C  eigenvector for RT1, giving the decomposition
C
C     [ CS1  SN1 ]  .  [ A B ]  .  [ CS1 -SN1 ]  =  [ RT1  0  ]
C     [-SN1  CS1 ]     [ B C ]     [ SN1  CS1 ]     [  0  RT2 ]
C
C  RT1 is accurate to a few ulps barring over/underflow.
C  RT2 may be inaccurate if there is massive cancellation in the
C     determinant A*C-B*B; higher precision or correctly rounded or
C     correctly truncated arithmetic would be needed to compute RT2
C     accurately in all cases.
C  CS1 and SN1 are accurate to a few ulps barring over/underflow.
C  Overflow is possible only if RT1 is within a factor of 5 of overflow.
C  Underflow is harmless if the input data is 0 or exceeds
C     underflow_threshold / macheps.
C
C  Arguments
C  =========
C
C  A      (input) REAL
C         The (1,1) entry of the input matrix.
C
C  B      (input) REAL
C         The (1,2) entry and the conjugate of the (2,1) entry of the
C         input matrix.
C
C  C      (input) REAL
C         The (2,2) entry of the input matrix.
C
C  RT1    (output) REAL
C         The eigenvalue of larger absolute value.
C
C  RT2    (output) REAL
C         The eigenvalue of smaller absolute value.
C
C  CS1    (output) REAL
C  SN1    (output) REAL
C         The vector (CS1, SN1) is a unit right eigenvector for RT1.
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
      DOUBLE PRECISION  A, B, C, CS1, RT1, RT2, SN1
C     .. Local Scalars ..
      DOUBLE PRECISION  AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     *                  TB, TN
      INTEGER           SGN1, SGN2
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
         SGN1 = -1
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE IF (SM.GT.ZERO) THEN
         RT1 = HALF*(SM+RT)
         SGN1 = 1
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
         SGN1 = 1
      END IF
C
C     Compute the eigenvector
C
      IF (DF.GE.ZERO) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS(CS)
      IF (ACS.GT.AB) THEN
         CT = -TB/CS
         SN1 = ONE/SQRT(ONE+CT*CT)
         CS1 = CT*SN1
      ELSE
         IF (AB.EQ.ZERO) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS/TB
            CS1 = ONE/SQRT(ONE+TN*TN)
            SN1 = TN*CS1
         END IF
      END IF
      IF (SGN1.EQ.SGN2) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
C
C     End of F07MDX (DLAEV2)
C
      END
