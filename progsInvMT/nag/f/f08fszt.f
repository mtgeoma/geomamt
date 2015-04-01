      SUBROUTINE F08FSZ(UPLO,N,A,LDA,D,E,TAU,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZHETD2(UPLO,N,A,LDA,D,E,TAU,INFO)
C
C  Purpose
C  =======
C
C  ZHETD2 reduces a complex Hermitian matrix A to real symmetric
C  tridiagonal form T by a unitary similarity transformation:
C  Q' * A * Q = T.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit, if UPLO = 'U', the diagonal and first superdiagonal
C          of A are overwritten by the corresponding elements of the
C          tridiagonal matrix T, and the elements above the first
C          superdiagonal, with the array TAU, represent the unitary
C          matrix Q as a product of elementary reflectors; if UPLO
C          = 'L', the diagonal and first subdiagonal of A are over-
C          written by the corresponding elements of the tridiagonal
C          matrix T, and the elements below the first subdiagonal, with
C          the array TAU, represent the unitary matrix Q as a product
C          of elementary reflectors. See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  D       (output) DOUBLE PRECISION array, dimension (N)
C          The diagonal elements of the tridiagonal matrix T:
C          D(i) = A(i,i).
C
C  E       (output) DOUBLE PRECISION array, dimension (N-1)
C          The off-diagonal elements of the tridiagonal matrix T:
C          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
C
C  TAU     (output) COMPLEX*16 array, dimension (N-1)
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(n-1) . . . H(2) H(1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
C  A(1:i-1,i+1), and tau in TAU(i).
C
C  If UPLO = 'L', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(1) H(2) . . . H(n-1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
C  and tau in TAU(i).
C
C  The contents of A on exit are illustrated by the following examples
C  with n = 5:
C
C  if UPLO = 'U':                       if UPLO = 'L':
C
C    (  d   e   v2  v3  v4 )              (  d                  )
C    (      d   e   v3  v4 )              (  e   d              )
C    (          d   e   v4 )              (  v1  e   d          )
C    (              d   e  )              (  v1  v2  e   d      )
C    (                  d  )              (  v1  v2  v3  e   d  )
C
C  where d and e denote diagonal and off-diagonal elements of T, and vi
C  denotes an element of the vector defining H(i).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO, HALF
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0,HALF=1.0D0/2.0D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*)
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA, TAUI
      INTEGER           I
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASV, ZAXPY, ZHEMV, ZHER2
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08FSZ/ZHETD2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
      IF (UPPER) THEN
C
C        Reduce the upper triangle of A
C
         A(N,N) = DBLE(A(N,N))
         DO 20 I = N - 1, 1, -1
C
C           Generate elementary reflector H(i) = I - tau * v * v'
C           to annihilate A(1:i-1,i+1)
C
            ALPHA = A(I,I+1)
            CALL F08ASV(I,ALPHA,A(1,I+1),1,TAUI)
            E(I) = ALPHA
C
            IF (TAUI.NE.ZERO) THEN
C
C              Apply H(i) from both sides to A(1:i,1:i)
C
               A(I,I+1) = ONE
C
C              Compute  x := tau * A * v  storing x in TAU(1:i)
C
               CALL ZHEMV(UPLO,I,TAUI,A,LDA,A(1,I+1),1,ZERO,TAU,1)
C
C              Compute  w := x - 1/2 * tau * (x'*v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC(I,TAU,1,A(1,I+1),1)
               CALL ZAXPY(I,ALPHA,A(1,I+1),1,TAU,1)
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w' - w * v'
C
               CALL ZHER2(UPLO,I,-ONE,A(1,I+1),1,TAU,1,A,LDA)
C
            END IF
            A(I,I+1) = E(I)
            D(I+1) = A(I+1,I+1)
            TAU(I) = TAUI
   20    CONTINUE
         D(1) = A(1,1)
      ELSE
C
C        Reduce the lower triangle of A
C
         A(1,1) = DBLE(A(1,1))
         DO 40 I = 1, N - 1
C
C           Generate elementary reflector H(i) = I - tau * v * v'
C           to annihilate A(i+2:n,i)
C
            ALPHA = A(I+1,I)
            CALL F08ASV(N-I,ALPHA,A(MIN(I+2,N),I),1,TAUI)
            E(I) = ALPHA
C
            IF (TAUI.NE.ZERO) THEN
C
C              Apply H(i) from both sides to A(i+1:n,i+1:n)
C
               A(I+1,I) = ONE
C
C              Compute  x := tau * A * v  storing y in TAU(i:n-1)
C
               CALL ZHEMV(UPLO,N-I,TAUI,A(I+1,I+1),LDA,A(I+1,I),1,ZERO,
     *                    TAU(I),1)
C
C              Compute  w := x - 1/2 * tau * (x'*v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC(N-I,TAU(I),1,A(I+1,I),1)
               CALL ZAXPY(N-I,ALPHA,A(I+1,I),1,TAU(I),1)
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w' - w * v'
C
               CALL ZHER2(UPLO,N-I,-ONE,A(I+1,I),1,TAU(I),1,A(I+1,I+1),
     *                    LDA)
C
            END IF
            A(I+1,I) = E(I)
            D(I) = A(I,I)
            TAU(I) = TAUI
   40    CONTINUE
         D(N) = A(N,N)
      END IF
C
      RETURN
C
C     End of F08FSZ (ZHETD2)
C
      END
