      SUBROUTINE F08CVY(M,N,A,LDA,TAU,WORK,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  ZGERQ2 computes an RQ factorization of a complex m by n matrix A:
C  A = R * Q.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the m by n matrix A.
C          On exit, if m <= n, the upper triangle of the subarray
C          A(1:m,n-m+1:n) contains the m by m upper triangular matrix R;
C          if m >= n, the elements on and above the (m-n)-th subdiagonal
C          contain the m by n upper trapezoidal matrix R; the remaining
C          elements, with the array TAU, represent the unitary matrix
C          Q as a product of elementary reflectors (see Further
C          Details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (M)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  Further Details
C  ===============
C
C  The matrix Q is represented as a product of elementary reflectors
C
C     Q = H(1)' H(2)' . . . H(k)', where k = min(m,n).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on
C  exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i).
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA
      INTEGER           I, K
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07FRY, F08ASV, F08ASW
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08CVY',-INFO)
         RETURN
      END IF
C
      K = MIN(M,N)
C
      DO 20 I = K, 1, -1
C
C        Generate elementary reflector H(i) to annihilate
C        A(m-k+i,1:n-k+i-1)
C
         CALL F07FRY(N-K+I,A(M-K+I,1),LDA)
         ALPHA = A(M-K+I,N-K+I)
         CALL F08ASV(N-K+I,ALPHA,A(M-K+I,1),LDA,TAU(I))
C
C        Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right
C
         A(M-K+I,N-K+I) = ONE
         CALL F08ASW('Right',M-K+I-1,N-K+I,A(M-K+I,1),LDA,TAU(I),A,LDA,
     *               WORK)
         A(M-K+I,N-K+I) = ALPHA
         CALL F07FRY(N-K+I-1,A(M-K+I,1),LDA)
   20 CONTINUE
      RETURN
C
C     End of F08CVY (ZGERQ2)
C
      END
