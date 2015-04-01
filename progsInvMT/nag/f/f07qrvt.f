      SUBROUTINE F07QRV(UPLO,N,ALPHA,X,INCX,AP)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZSPR(UPLO,N,ALPHA,X,INCX,AP)
C
C  Purpose
C  =======
C
C  ZSPR    performs the symmetric rank 1 operation
C
C     A := alpha*x*conjg( x' ) + A,
C
C  where alpha is a real scalar, x is an n element vector and A is an
C  n by n symmetric matrix, supplied in packed form.
C
C  Arguments
C  ==========
C
C  UPLO   - CHARACTER*1
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the matrix A is supplied in the packed
C           array AP as follows:
C
C              UPLO = 'U' or 'u'   The upper triangular part of A is
C                                  supplied in AP.
C
C              UPLO = 'L' or 'l'   The lower triangular part of A is
C                                  supplied in AP.
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
C  AP     - COMPLEX array, dimension at least
C           ( ( N*( N + 1 ) )/2 ).
C           Before entry with  UPLO = 'U' or 'u', the array AP must
C           contain the upper triangular part of the symmetric matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
C           and a( 2, 2 ) respectively, and so on. On exit, the array
C           AP is overwritten by the upper triangular part of the
C           updated matrix.
C           Before entry with UPLO = 'L' or 'l', the array AP must
C           contain the lower triangular part of the symmetric matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
C           and a( 3, 1 ) respectively, and so on. On exit, the array
C           AP is overwritten by the lower triangular part of the
C           updated matrix.
C           Note that the imaginary parts of the diagonal elements need
C           not be set, they are assumed to be zero, and on exit they
C           are set to zero.
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
      INTEGER           INCX, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), X(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, IX, J, JX, K, KK, KX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
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
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07QRV/ZSPR',INFO)
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
C     Start the operations. In this version the elements of the array AP
C     are accessed sequentially with one pass through AP.
C
      KK = 1
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Form  A  when upper triangle is stored in AP.
C
         IF (INCX.EQ.1) THEN
            DO 40 J = 1, N
               IF (X(J).NE.ZERO) THEN
                  TEMP = ALPHA*X(J)
                  K = KK
                  DO 20 I = 1, J - 1
                     AP(K) = AP(K) + X(I)*TEMP
                     K = K + 1
   20             CONTINUE
                  AP(KK+J-1) = AP(KK+J-1) + X(J)*TEMP
               ELSE
                  AP(KK+J-1) = AP(KK+J-1)
               END IF
               KK = KK + J
   40       CONTINUE
         ELSE
            JX = KX
            DO 80 J = 1, N
               IF (X(JX).NE.ZERO) THEN
                  TEMP = ALPHA*X(JX)
                  IX = KX
                  DO 60 K = KK, KK + J - 2
                     AP(K) = AP(K) + X(IX)*TEMP
                     IX = IX + INCX
   60             CONTINUE
                  AP(KK+J-1) = AP(KK+J-1) + X(JX)*TEMP
               ELSE
                  AP(KK+J-1) = AP(KK+J-1)
               END IF
               JX = JX + INCX
               KK = KK + J
   80       CONTINUE
         END IF
      ELSE
C
C        Form  A  when lower triangle is stored in AP.
C
         IF (INCX.EQ.1) THEN
            DO 120 J = 1, N
               IF (X(J).NE.ZERO) THEN
                  TEMP = ALPHA*X(J)
                  AP(KK) = AP(KK) + TEMP*X(J)
                  K = KK + 1
                  DO 100 I = J + 1, N
                     AP(K) = AP(K) + X(I)*TEMP
                     K = K + 1
  100             CONTINUE
               ELSE
                  AP(KK) = AP(KK)
               END IF
               KK = KK + N - J + 1
  120       CONTINUE
         ELSE
            JX = KX
            DO 160 J = 1, N
               IF (X(JX).NE.ZERO) THEN
                  TEMP = ALPHA*X(JX)
                  AP(KK) = AP(KK) + TEMP*X(JX)
                  IX = JX
                  DO 140 K = KK + 1, KK + N - J
                     IX = IX + INCX
                     AP(K) = AP(K) + X(IX)*TEMP
  140             CONTINUE
               ELSE
                  AP(KK) = AP(KK)
               END IF
               JX = JX + INCX
               KK = KK + N - J + 1
  160       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F07QRV (ZSPR)
C
      END
