      SUBROUTINE F07NRV(UPLO,N,ALPHA,X,INCX,A,LDA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
C
C  Purpose
C  =======
C
C  ZSYR   performs the symmetric rank 1 operation
C
C     A := alpha*x*( x' ) + A,
C
C  where alpha is a real scalar, x is an n element vector and A is an
C  n by n symmetric matrix.
C
C  Arguments
C  ==========
C
C  UPLO   - CHARACTER*1
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - COMPLEX array, dimension at least
C           ( 1 + ( N - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the N-
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  A      - COMPLEX array, dimension( LDA, N )
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, N ).
C           Unchanged on exit.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Rice University
C
C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA
      INTEGER           INCX, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), X(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, IX, J, JX, KX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF ( .NOT. (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
     *    .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = 1
      ELSE IF (N.LT.0) THEN
         INFO = 2
      ELSE IF (INCX.EQ.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = 7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07NRV/ZSYR',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
C
C     Set the start point in X if the increment is not unity.
C
      IF (INCX.LE.0) THEN
         KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Form  A  when A is stored in upper triangle.
C
         IF (INCX.EQ.1) THEN
            DO 40 J = 1, N
               IF (X(J).NE.ZERO) THEN
                  TEMP = ALPHA*X(J)
                  DO 20 I = 1, J
                     A(I,J) = A(I,J) + X(I)*TEMP
   20             CONTINUE
               END IF
   40       CONTINUE
         ELSE
            JX = KX
            DO 80 J = 1, N
               IF (X(JX).NE.ZERO) THEN
                  TEMP = ALPHA*X(JX)
                  IX = KX
                  DO 60 I = 1, J
                     A(I,J) = A(I,J) + X(IX)*TEMP
                     IX = IX + INCX
   60             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in lower triangle.
C
         IF (INCX.EQ.1) THEN
            DO 120 J = 1, N
               IF (X(J).NE.ZERO) THEN
                  TEMP = ALPHA*X(J)
                  DO 100 I = J, N
                     A(I,J) = A(I,J) + X(I)*TEMP
  100             CONTINUE
               END IF
  120       CONTINUE
         ELSE
            JX = KX
            DO 160 J = 1, N
               IF (X(JX).NE.ZERO) THEN
                  TEMP = ALPHA*X(JX)
                  IX = JX
                  DO 140 I = J, N
                     A(I,J) = A(I,J) + X(IX)*TEMP
                     IX = IX + INCX
  140             CONTINUE
               END IF
               JX = JX + INCX
  160       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F07NRV (ZSYR)
C
      END
