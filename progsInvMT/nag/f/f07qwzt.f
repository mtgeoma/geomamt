      SUBROUTINE F07QWZ(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
C
C  Purpose
C  =======
C
C  ZSPMV  performs the matrix-vector operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n hermitian matrix, supplied in packed form.
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
C  AP     - COMPLEX array, dimension at least
C           ( ( N*( N + 1 ) )/2 ).
C           Before entry with UPLO = 'U' or 'u', the array AP must
C           contain the upper triangular part of the hermitian matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
C           and a( 2, 2 ) respectively, and so on.
C           Before entry with UPLO = 'L' or 'l', the array AP must
C           contain the lower triangular part of the hermitian matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
C           and a( 3, 1 ) respectively, and so on.
C           Note that the imaginary parts of the diagonal elements need
C           not be set and are assumed to be zero.
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
C  BETA   - COMPLEX
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - COMPLEX array, dimension at least
C           ( 1 + ( N - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Rice University
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           INCX, INCY, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), X(*), Y(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP1, TEMP2
      INTEGER           I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
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
         INFO = 6
      ELSE IF (INCY.EQ.0) THEN
         INFO = 9
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07QWZ/ZSPMV',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO) .AND. (BETA.EQ.ONE))) RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF (INCX.GT.0) THEN
         KX = 1
      ELSE
         KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
         KY = 1
      ELSE
         KY = 1 - (N-1)*INCY
      END IF
C
C     Start the operations. In this version the elements of the array AP
C     are accessed sequentially with one pass through AP.
C
C     First form  y := beta*y.
C
      IF (BETA.NE.ONE) THEN
         IF (INCY.EQ.1) THEN
            IF (BETA.EQ.ZERO) THEN
               DO 20 I = 1, N
                  Y(I) = ZERO
   20          CONTINUE
            ELSE
               DO 40 I = 1, N
                  Y(I) = BETA*Y(I)
   40          CONTINUE
            END IF
         ELSE
            IY = KY
            IF (BETA.EQ.ZERO) THEN
               DO 60 I = 1, N
                  Y(IY) = ZERO
                  IY = IY + INCY
   60          CONTINUE
            ELSE
               DO 80 I = 1, N
                  Y(IY) = BETA*Y(IY)
                  IY = IY + INCY
   80          CONTINUE
            END IF
         END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      KK = 1
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Form  y  when AP contains the upper triangle.
C
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
            DO 120 J = 1, N
               TEMP1 = ALPHA*X(J)
               TEMP2 = ZERO
               K = KK
               DO 100 I = 1, J - 1
                  Y(I) = Y(I) + TEMP1*AP(K)
                  TEMP2 = TEMP2 + AP(K)*X(I)
                  K = K + 1
  100          CONTINUE
               Y(J) = Y(J) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
               KK = KK + J
  120       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 160 J = 1, N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = ZERO
               IX = KX
               IY = KY
               DO 140 K = KK, KK + J - 2
                  Y(IY) = Y(IY) + TEMP1*AP(K)
                  TEMP2 = TEMP2 + AP(K)*X(IX)
                  IX = IX + INCX
                  IY = IY + INCY
  140          CONTINUE
               Y(JY) = Y(JY) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + J
  160       CONTINUE
         END IF
      ELSE
C
C        Form  y  when AP contains the lower triangle.
C
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
            DO 200 J = 1, N
               TEMP1 = ALPHA*X(J)
               TEMP2 = ZERO
               Y(J) = Y(J) + TEMP1*AP(KK)
               K = KK + 1
               DO 180 I = J + 1, N
                  Y(I) = Y(I) + TEMP1*AP(K)
                  TEMP2 = TEMP2 + AP(K)*X(I)
                  K = K + 1
  180          CONTINUE
               Y(J) = Y(J) + ALPHA*TEMP2
               KK = KK + (N-J+1)
  200       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 240 J = 1, N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = ZERO
               Y(JY) = Y(JY) + TEMP1*AP(KK)
               IX = JX
               IY = JY
               DO 220 K = KK + 1, KK + N - J
                  IX = IX + INCX
                  IY = IY + INCY
                  Y(IY) = Y(IY) + TEMP1*AP(K)
                  TEMP2 = TEMP2 + AP(K)*X(IX)
  220          CONTINUE
               Y(JY) = Y(JY) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + (N-J+1)
  240       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F07QWZ (ZSPMV)
C
      END
