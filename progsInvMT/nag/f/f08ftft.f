      SUBROUTINE F08FTF(UPLO,N,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZUNGTR(UPLO,N,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNGTR generates a complex unitary matrix Q which is defined as the
C  product of n-1 elementary reflectors of order n, as returned by
C  ZHETRD:
C
C  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
C
C  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangle of the array A
C          holds details of the elementary reflectors, as returned by
C          ZHETRD:
C          = 'U': Upper triangle;
C          = 'L': Lower triangle.
C
C  N       (input) INTEGER
C          The order of the matrix Q. N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the vectors which define the elementary reflectors,
C          as returned by ZHETRD.
C          On exit, the n by n unitary matrix Q.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= N.
C
C  TAU     (input) COMPLEX*16 array, dimension (N-1)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZHETRD.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= N-1.
C          For optimum performance LWORK should be at least (N-1)*NB,
C          where NB is the optimal blocksize.
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
      INTEGER           INFO, LDA, LWORK, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IINFO, J
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08CTY, ZUNGQR
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
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.MAX(1,N-1)) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08FTF/ZUNGTR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
      IF (UPPER) THEN
C
C        Q was determined by a call to ZHETRD with UPLO = 'U'
C
C        Shift the vectors which define the elementary reflectors one
C        column to the left, and set the last row and column of Q to
C        those of the unit matrix
C
         DO 40 J = 1, N - 1
            DO 20 I = 1, J - 1
               A(I,J) = A(I,J+1)
   20       CONTINUE
            A(N,J) = ZERO
   40    CONTINUE
         DO 60 I = 1, N - 1
            A(I,N) = ZERO
   60    CONTINUE
         A(N,N) = ONE
C
C        Generate Q(1:n-1,1:n-1)
C
         CALL F08CTY(N-1,N-1,N-1,A,LDA,TAU,WORK,LWORK,IINFO)
C
      ELSE
C
C        Q was determined by a call to ZHETRD with UPLO = 'L'.
C
C        Shift the vectors which define the elementary reflectors one
C        column to the right, and set the first row and column of Q to
C        those of the unit matrix
C
         DO 100 J = N, 2, -1
            A(1,J) = ZERO
            DO 80 I = J + 1, N
               A(I,J) = A(I,J-1)
   80       CONTINUE
  100    CONTINUE
         A(1,1) = ONE
         DO 120 I = 2, N
            A(I,1) = ZERO
  120    CONTINUE
         IF (N.GT.1) THEN
C
C           Generate Q(2:n,2:n)
C
            CALL ZUNGQR(N-1,N-1,N-1,A(2,2),LDA,TAU,WORK,LWORK,IINFO)
         ELSE
            WORK(1) = 1
         END IF
      END IF
      RETURN
C
C     End of F08FTF (ZUNGTR)
C
      END
