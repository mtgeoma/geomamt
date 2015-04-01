      SUBROUTINE F07MRX(A,B,C,RT1,RT2,CS1,SN1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLAEV2(A,B,C,RT1,RT2,CS1,SN1)
C
C  Purpose
C  =======
C
C  ZLAEV2 computes the eigendecomposition of a 2 by 2 Hermitian matrix:
C     ((A,B);(CONJG(B),C))
C
C  RT1 is the eigenvalue of larger absolute value, RT2 is the eigenvalue
C  of smaller absolute value, and (CS1,SN1) is the unit right
C  eigenvector for RT1, giving the decomposition
C
C  [ CS1  CONJG(SN1) ].[    A     B ].[ CS1 -CONJG(SN1) ] =  [ RT1  0  ]
C  [-SN1     CS1     ] [ CONJG(B) C ] [ SN1     CS1     ]    [  0  RT2 ]
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
C  A      (input) COMPLEX
C         The (1,1) entry of the input matrix.
C
C  B      (input) COMPLEX
C         The (1,2) entry and the conjugate of the (2,1) entry of the
C         input matrix.
C
C  C      (input) COMPLEX
C         The (2,2) entry of the input matrix.
C
C  RT1    (output) REAL
C         The eigenvalue of larger absolute value.
C
C  RT2    (output) REAL
C         The eigenvalue of smaller absolute value.
C
C  CS1    (output) REAL
C  SN1    (output) COMPLEX
C         The vector (CS1, SN1) is a unit right eigenvector for RT1.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
C     .. Scalar Arguments ..
      COMPLEX*16        A, B, C, SN1
      DOUBLE PRECISION  CS1, RT1, RT2
C     .. Local Scalars ..
      COMPLEX*16        W
      DOUBLE PRECISION  T
C     .. External Subroutines ..
      EXTERNAL          F07MDX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCONJG, DBLE
C     .. Executable Statements ..
C
      IF (ABS(B).EQ.ZERO) THEN
         W = ONE
      ELSE
         W = DCONJG(B)/ABS(B)
      END IF
      CALL F07MDX(DBLE(A),ABS(B),DBLE(C),RT1,RT2,CS1,T)
      SN1 = W*T
      RETURN
C
C     End of F07MRX (ZLAEV2)
C
      END
