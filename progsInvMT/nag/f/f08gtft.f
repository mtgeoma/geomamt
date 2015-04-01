      SUBROUTINE F08GTF(UPLO,N,AP,TAU,Q,LDQ,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZUPGTR(UPLO,N,AP,TAU,Q,LDQ,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZUPGTR generates a complex unitary matrix Q which is defined as the
C  product of n-1 elementary reflectors of order n, as returned by
C  ZHPTRD using packed storage:
C
C  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
C
C  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies the storage scheme used in the previous call of
C          ZHPTRD:
C          = 'U': Upper triangular packed storage;
C          = 'L': Lower triangular packed storage.
C
C  N       (input) INTEGER
C          The order of the matrix Q. N >= 0.
C
C  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
C          The vectors which define the elementary reflectors, as
C          returned by ZHPTRD.
C
C  TAU     (input) COMPLEX*16 array, dimension (N-1)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZHPTRD.
C
C  Q       (output) COMPLEX*16 array, dimension (LDQ,N)
C          The n-by-n unitary matrix Q.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q. LDQ >= max(1,N).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N-1)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDQ, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), Q(LDQ,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      INTEGER           I, IINFO, IJ, J
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ATZ, F08CTZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDQ.LT.MAX(1,N)) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08GTF/ZUPGTR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Q was determined by a call to ZHPTRD with UPLO = 'U'
C
C        Unpack the vectors which define the elementary reflectors and
C        set the last row and column of Q equal to those of the unit
C        matrix
C
         IJ = 2
         DO 40 J = 1, N - 1
            DO 20 I = 1, J - 1
               Q(I,J) = AP(IJ)
               IJ = IJ + 1
   20       CONTINUE
            IJ = IJ + 2
            Q(N,J) = ZERO
   40    CONTINUE
         DO 60 I = 1, N - 1
            Q(I,N) = ZERO
   60    CONTINUE
         Q(N,N) = ONE
C
C        Generate Q(1:n-1,1:n-1)
C
         CALL F08CTZ(N-1,N-1,N-1,Q,LDQ,TAU,WORK,IINFO)
C
      ELSE
C
C        Q was determined by a call to ZHPTRD with UPLO = 'L'.
C
C        Unpack the vectors which define the elementary reflectors and
C        set the first row and column of Q equal to those of the unit
C        matrix
C
         Q(1,1) = ONE
         DO 80 I = 2, N
            Q(I,1) = ZERO
   80    CONTINUE
         IJ = 3
         DO 120 J = 2, N
            Q(1,J) = ZERO
            DO 100 I = J + 1, N
               Q(I,J) = AP(IJ)
               IJ = IJ + 1
  100       CONTINUE
            IJ = IJ + 2
  120    CONTINUE
         IF (N.GT.1) THEN
C
C           Generate Q(2:n,2:n)
C
            CALL F08ATZ(N-1,N-1,N-1,Q(2,2),LDQ,TAU,WORK,IINFO)
         END IF
      END IF
      RETURN
C
C     End of F08GTF (ZUPGTR)
C
      END
