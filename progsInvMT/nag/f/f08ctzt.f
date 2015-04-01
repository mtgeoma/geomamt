      SUBROUTINE F08CTZ(M,N,K,A,LDA,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZUNG2L(M,N,K,A,LDA,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNG2L generates an m by n complex matrix Q with orthonormal columns,
C  which is defined as the last n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(k) . . . H(2) H(1)
C
C  as returned by ZGEQLF.
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
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the (n-k+i)-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by ZGEQLF in the last k columns of its array
C          argument A.
C          On exit, the m-by-n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEQLF.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N)
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
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      INTEGER           I, II, J, L
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASW, ZSCAL
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
         CALL F06AAZ('F08CTZ/ZUNG2L',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
C     Initialise columns 1:n-k to columns of the unit matrix
C
      DO 40 J = 1, N - K
         DO 20 L = 1, M
            A(L,J) = ZERO
   20    CONTINUE
         A(M-N+J,J) = ONE
   40 CONTINUE
C
      DO 80 I = 1, K
         II = N - K + I
C
C        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
C
         A(M-N+II,II) = ONE
         CALL F08ASW('Left',M-N+II,II-1,A(1,II),1,TAU(I),A,LDA,WORK)
         CALL ZSCAL(M-N+II-1,-TAU(I),A(1,II),1)
         A(M-N+II,II) = ONE - TAU(I)
C
C        Set A(m-k+i+1:m,n-k+i) to zero
C
         DO 60 L = M - N + II + 1, M
            A(L,II) = ZERO
   60    CONTINUE
   80 CONTINUE
      RETURN
C
C     End of F08CTZ (ZUNG2L)
C
      END
