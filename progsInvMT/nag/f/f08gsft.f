      SUBROUTINE F08GSF(UPLO,N,AP,D,E,TAU,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZHPTRD(UPLO,N,AP,D,E,TAU,INFO)
C
C  Purpose
C  =======
C
C  ZHPTRD reduces a complex Hermitian matrix A stored in packed form to
C  real symmetric tridiagonal form T by a unitary similarity
C  transformation: Q' * A * Q = T.
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
C  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
C          On entry, the upper or lower triangle of the Hermitian matrix
C          A, packed columnwise in a linear array.  The j-th column of A
C          is stored in the array AP as follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
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
C          = 0:  successful exit.
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
C  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
C  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
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
C  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
C  overwriting A(i+2:n,i), and tau is stored in TAU(i).
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
      INTEGER           INFO, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), TAU(*)
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA, TAUI
      INTEGER           I, I1, I1I1, II
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASV, ZAXPY, ZHPMV, ZHPR2
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
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
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08GSF/ZHPTRD',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
      IF (UPPER) THEN
C
C        Reduce the upper triangle of A.
C        I1 is the index in AP of A(1,I+1).
C
         I1 = N*(N-1)/2 + 1
         AP(I1+N-1) = DBLE(AP(I1+N-1))
         DO 20 I = N - 1, 1, -1
C
C           Generate elementary reflector H(i) = I - tau * v * v'
C           to annihilate A(1:i-1,i+1)
C
            ALPHA = AP(I1+I-1)
            CALL F08ASV(I,ALPHA,AP(I1),1,TAUI)
            E(I) = ALPHA
C
            IF (TAUI.NE.ZERO) THEN
C
C              Apply H(i) from both sides to A(1:i,1:i)
C
               AP(I1+I-1) = ONE
C
C              Compute  y := tau * A * v  storing y in TAU(1:i)
C
               CALL ZHPMV(UPLO,I,TAUI,AP,AP(I1),1,ZERO,TAU,1)
C
C              Compute  w := y - 1/2 * tau * (y'*v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC(I,TAU,1,AP(I1),1)
               CALL ZAXPY(I,ALPHA,AP(I1),1,TAU,1)
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w' - w * v'
C
               CALL ZHPR2(UPLO,I,-ONE,AP(I1),1,TAU,1,AP)
C
            END IF
            AP(I1+I-1) = E(I)
            D(I+1) = AP(I1+I)
            TAU(I) = TAUI
            I1 = I1 - I
   20    CONTINUE
         D(1) = AP(1)
      ELSE
C
C        Reduce the lower triangle of A. II is the index in AP of
C        A(i,i) and I1I1 is the index of A(i+1,i+1).
C
         II = 1
         AP(1) = DBLE(AP(1))
         DO 40 I = 1, N - 1
            I1I1 = II + N - I + 1
C
C           Generate elementary reflector H(i) = I - tau * v * v'
C           to annihilate A(i+2:n,i)
C
            ALPHA = AP(II+1)
            CALL F08ASV(N-I,ALPHA,AP(II+2),1,TAUI)
            E(I) = ALPHA
C
            IF (TAUI.NE.ZERO) THEN
C
C              Apply H(i) from both sides to A(i+1:n,i+1:n)
C
               AP(II+1) = ONE
C
C              Compute  y := tau * A * v  storing y in TAU(i:n-1)
C
               CALL ZHPMV(UPLO,N-I,TAUI,AP(I1I1),AP(II+1),1,ZERO,TAU(I),
     *                    1)
C
C              Compute  w := y - 1/2 * tau * (y'*v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC(N-I,TAU(I),1,AP(II+1),1)
               CALL ZAXPY(N-I,ALPHA,AP(II+1),1,TAU(I),1)
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w' - w * v'
C
               CALL ZHPR2(UPLO,N-I,-ONE,AP(II+1),1,TAU(I),1,AP(I1I1))
C
            END IF
            AP(II+1) = E(I)
            D(I) = AP(II)
            TAU(I) = TAUI
            II = I1I1
   40    CONTINUE
         D(N) = AP(II)
      END IF
C
      RETURN
C
C     End of F08GSF (ZHPTRD)
C
      END
