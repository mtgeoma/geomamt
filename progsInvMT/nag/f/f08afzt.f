      SUBROUTINE F08AFZ(M,N,K,A,LDA,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DORG2R(M,N,K,A,LDA,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  DORG2R generates an m by n real matrix Q with orthonormal columns,
C  which is defined as the first n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(1) H(2) . . . H(k)
C
C  as returned by DGEQRF.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q. M >= N >= 0.
C
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines the
C          matrix Q. N >= K >= 0.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the i-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by DGEQRF in the first k columns of its array
C          argument A.
C          On exit, the m-by-n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) DOUBLE PRECISION array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by DGEQRF.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument has an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      INTEGER           I, J, L
C     .. External Subroutines ..
      EXTERNAL          DSCAL, F06AAZ, F08AEW
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0 .OR. N.GT.M) THEN
         INFO = -2
      ELSE IF (K.LT.0 .OR. K.GT.N) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08AFZ/DORG2R',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
C     Initialise columns k+1:n to columns of the unit matrix
C
      DO 40 J = K + 1, N
         DO 20 L = 1, M
            A(L,J) = ZERO
   20    CONTINUE
         A(J,J) = ONE
   40 CONTINUE
C
      DO 80 I = K, 1, -1
C
C        Apply H(i) to A(i:m,i:n) from the left
C
         IF (I.LT.N) THEN
            A(I,I) = ONE
            CALL F08AEW('Left',M-I+1,N-I,A(I,I),1,TAU(I),A(I,I+1),LDA,
     *                  WORK)
         END IF
         IF (I.LT.M) CALL DSCAL(M-I,-TAU(I),A(I+1,I),1)
         A(I,I) = ONE - TAU(I)
C
C        Set A(1:i-1,i) to zero
C
         DO 60 L = 1, I - 1
            A(L,I) = ZERO
   60    CONTINUE
   80 CONTINUE
      RETURN
C
C     End of F08AFZ (DORG2R)
C
      END
