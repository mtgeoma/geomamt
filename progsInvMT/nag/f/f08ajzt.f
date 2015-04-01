      SUBROUTINE F08AJZ(M,N,K,A,LDA,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DORGL2(M,N,K,A,LDA,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  DORGL2 generates an m by n real matrix Q with orthonormal rows,
C  which is defined as the first m rows of a product of k elementary
C  reflectors of order n
C
C        Q  =  H(k) . . . H(2) H(1)
C
C  as returned by DGELQF.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q. N >= M.
C
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines the
C          matrix Q. M >= K >= 0.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the i-th row must contain the vector which defines
C          the elementary reflector H(i), for i = 1,2,...,k, as returned
C          by DGELQF in the first k rows of its array argument A.
C          On exit, the m-by-n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) DOUBLE PRECISION array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by DGELQF.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
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
      ELSE IF (N.LT.M) THEN
         INFO = -2
      ELSE IF (K.LT.0 .OR. K.GT.M) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08AJZ/DORGL2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.LE.0) RETURN
C
      IF (K.LT.M) THEN
C
C        Initialise rows k+1:m to rows of the unit matrix
C
         DO 40 J = 1, N
            DO 20 L = K + 1, M
               A(L,J) = ZERO
   20       CONTINUE
            IF (J.GT.K .AND. J.LE.M) A(J,J) = ONE
   40    CONTINUE
      END IF
C
      DO 80 I = K, 1, -1
C
C        Apply H(i) to A(i:m,i:n) from the right
C
         IF (I.LT.N) THEN
            IF (I.LT.M) THEN
               A(I,I) = ONE
               CALL F08AEW('Right',M-I,N-I+1,A(I,I),LDA,TAU(I),A(I+1,I),
     *                     LDA,WORK)
            END IF
            CALL DSCAL(N-I,-TAU(I),A(I,I+1),LDA)
         END IF
         A(I,I) = ONE - TAU(I)
C
C        Set A(1:i-1,i) to zero
C
         DO 60 L = 1, I - 1
            A(I,L) = ZERO
   60    CONTINUE
   80 CONTINUE
      RETURN
C
C     End of F08AJZ (DORGL2)
C
      END
