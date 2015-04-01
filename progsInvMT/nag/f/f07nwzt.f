      SUBROUTINE F07NWZ(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
C
C  Purpose
C  =======
C
C  ZSYMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric matrix.
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
C  A      - COMPLEX array, dimension( LDA, N )
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced.
C           Note that the imaginary parts of the diagonal elements need
C           not be set and are assumed to be zero.
C           Unchanged on exit.
C
C  LDA    - INTEGER
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, N ).
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
      INTEGER           INCX, INCY, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), X(*), Y(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP1, TEMP2
      INTEGER           I, INFO, IX, IY, J, JX, JY, KX, KY
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
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = 5
      ELSE IF (INCX.EQ.0) THEN
         INFO = 7
      ELSE IF (INCY.EQ.0) THEN
         INFO = 10
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07NWZ/ZSYMV',INFO)
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
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
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
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Form  y  when A is stored in upper triangle.
C
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
            DO 120 J = 1, N
               TEMP1 = ALPHA*X(J)
               TEMP2 = ZERO
               DO 100 I = 1, J - 1
                  Y(I) = Y(I) + TEMP1*A(I,J)
                  TEMP2 = TEMP2 + A(I,J)*X(I)
  100          CONTINUE
               Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
  120       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 160 J = 1, N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = ZERO
               IX = KX
               IY = KY
               DO 140 I = 1, J - 1
                  Y(IY) = Y(IY) + TEMP1*A(I,J)
                  TEMP2 = TEMP2 + A(I,J)*X(IX)
                  IX = IX + INCX
                  IY = IY + INCY
  140          CONTINUE
               Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
  160       CONTINUE
         END IF
      ELSE
C
C        Form  y  when A is stored in lower triangle.
C
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
            DO 200 J = 1, N
               TEMP1 = ALPHA*X(J)
               TEMP2 = ZERO
               Y(J) = Y(J) + TEMP1*A(J,J)
               DO 180 I = J + 1, N
                  Y(I) = Y(I) + TEMP1*A(I,J)
                  TEMP2 = TEMP2 + A(I,J)*X(I)
  180          CONTINUE
               Y(J) = Y(J) + ALPHA*TEMP2
  200       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 240 J = 1, N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = ZERO
               Y(JY) = Y(JY) + TEMP1*A(J,J)
               IX = JX
               IY = JY
               DO 220 I = J + 1, N
                  IX = IX + INCX
                  IY = IY + INCY
                  Y(IY) = Y(IY) + TEMP1*A(I,J)
                  TEMP2 = TEMP2 + A(I,J)*X(IX)
  220          CONTINUE
               Y(JY) = Y(JY) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
  240       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F07NWZ (ZSYMV)
C
      END
