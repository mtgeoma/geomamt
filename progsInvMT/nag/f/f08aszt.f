      SUBROUTINE F08ASZ(M,N,A,LDA,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZGEQR2(M,N,A,LDA,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZGEQR2 computes a QR factorization of a complex m by n matrix A:
C  A = Q * R.
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
C          On exit, the elements on and above the diagonal of the array
C          contain the min(m,n) by n upper trapezoidal matrix R (R is
C          upper triangular if m >= n); the elements below the diagonal,
C          with the array TAU, represent the unitary matrix Q as a
C          product of elementary reflectors (see Further Details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N)
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
C     Q = H(1) H(2) . . . H(k), where k = min(m,n).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
C  and tau in TAU(i).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA
      INTEGER           I, K
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASV, F08ASW
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX, MIN
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
         CALL F06AAZ('F08ASZ/ZGEQR2',-INFO)
         RETURN
      END IF
C
      K = MIN(M,N)
C
      DO 20 I = 1, K
C
C        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
C
         CALL F08ASV(M-I+1,A(I,I),A(MIN(I+1,M),I),1,TAU(I))
         IF (I.LT.N) THEN
C
C           Apply H(i)' to A(i:m,i+1:n) from the left
C
            ALPHA = A(I,I)
            A(I,I) = ONE
            CALL F08ASW('Left',M-I+1,N-I,A(I,I),1,DCONJG(TAU(I)),
     *                  A(I,I+1),LDA,WORK)
            A(I,I) = ALPHA
         END IF
   20 CONTINUE
      RETURN
C
C     End of F08ASZ (ZGEQR2)
C
      END
