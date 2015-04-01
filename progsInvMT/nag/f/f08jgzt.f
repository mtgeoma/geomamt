      SUBROUTINE F08JGZ(N,D,E,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DPTTRF(N,D,E,INFO)
C
C  Purpose
C  =======
C
C  DPTTRF computes the L*D*L' factorization of a real symmetric
C  positive definite tridiagonal matrix A.  The factorization may also
C  be regarded as having the form A = U'*D*U.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, the n diagonal elements of the tridiagonal matrix
C          A.  On exit, the n diagonal elements of the diagonal matrix
C          D from the L*D*L' factorization of A.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          On entry, the (n-1) subdiagonal elements of the tridiagonal
C          matrix A.  On exit, the (n-1) subdiagonal elements of the
C          unit bidiagonal factor L from the L*D*L' factorization of A.
C          E can also be regarded as the superdiagonal of the unit
C          bidiagonal factor U from the U'*D*U factorization of A.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the leading minor of order k is not
C               positive definite; if k < N, the factorization could not
C               be completed, while if k = N, the factorization was
C               completed, but D(N) = 0.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DI, EI
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF (N.LT.0) THEN
         INFO = -1
         CALL F06AAZ('F08JGZ/DPTTRF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Compute the L*D*L' (or U'*D*U) factorization of A.
C
      DO 20 I = 1, N - 1
C
C        Drop out of the loop if d(i) <= 0: the matrix is not positive
C        definite.
C
         DI = D(I)
         IF (DI.LE.ZERO) GO TO 40
C
C        Solve for e(i) and d(i+1).
C
         EI = E(I)
         E(I) = EI/DI
         D(I+1) = D(I+1) - E(I)*EI
   20 CONTINUE
C
C     Check d(n) for positive definiteness.
C
      I = N
      IF (D(I).GT.ZERO) GO TO 60
C
   40 CONTINUE
      INFO = I
C
   60 CONTINUE
      RETURN
C
C     End of F08JGZ (DPTTRF)
C
      END
