      SUBROUTINE F06ZFF(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZTRMM  performs one of the matrix-matrix operations
C
C     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
C
C  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
C  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C
C     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C           On entry,  SIDE specifies whether  op( A ) multiplies B from
C           the left or right as follows:
C
C              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
C
C              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
C
C           Unchanged on exit.
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix A is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n'   op( A ) = A.
C
C              TRANSA = 'T' or 't'   op( A ) = A'.
C
C              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit triangular
C           as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of B. M must be at
C           least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of B.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
C           zero then  A is not referenced and  B need not be set before
C           entry.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
C           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
C           upper triangular part of the array  A must contain the upper
C           triangular matrix  and the strictly lower triangular part of
C           A is not referenced.
C           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
C           lower triangular part of the array  A must contain the lower
C           triangular matrix  and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C           A  are not referenced either,  but are assumed to be  unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C           then LDA must be at least max( 1, n ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, n ).
C           Before entry,  the leading  m by n part of the array  B must
C           contain the matrix  B,  and  on exit  is overwritten  by the
C           transformed matrix.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   LDB  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
      ENTRY             ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,
     *                  LDB)
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA
      INTEGER           LDA, LDB, M, N
      CHARACTER*1       DIAG, SIDE, TRANSA, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, J, K, NROWA
      LOGICAL           LSIDE, NOCONJ, NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      LSIDE = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t')
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
      INFO = 0
      IF (( .NOT. LSIDE) .AND. ( .NOT. (SIDE.EQ.'R' .OR. SIDE.EQ.'r')))
     *    THEN
         INFO = 1
      ELSE IF (( .NOT. UPPER) .AND. ( .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.
     *         'l'))) THEN
         INFO = 2
      ELSE IF (( .NOT. (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n'))
     *         .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))
     *         .AND. ( .NOT. (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c'))) THEN
         INFO = 3
      ELSE IF (( .NOT. (DIAG.EQ.'U' .OR. DIAG.EQ.'u'))
     *         .AND. ( .NOT. (DIAG.EQ.'N' .OR. DIAG.EQ.'n'))) THEN
         INFO = 4
      ELSE IF (M.LT.0) THEN
         INFO = 5
      ELSE IF (N.LT.0) THEN
         INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZFF/ZTRMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF (N.EQ.0) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         DO 40 J = 1, N
            DO 20 I = 1, M
               B(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
         RETURN
      END IF
C
C     Start the operations.
C
      IF (LSIDE) THEN
         IF ((TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')) THEN
C
C           Form  B := alpha*A*B.
C
            IF (UPPER) THEN
               DO 100 J = 1, N
                  DO 80 K = 1, M
                     IF (B(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*B(K,J)
                        DO 60 I = 1, K - 1
                           B(I,J) = B(I,J) + TEMP*A(I,K)
   60                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP*A(K,K)
                        B(K,J) = TEMP
                     END IF
   80             CONTINUE
  100          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 140 K = M, 1, -1
                     IF (B(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*B(K,J)
                        B(K,J) = TEMP
                        IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                        DO 120 I = K + 1, M
                           B(I,J) = B(I,J) + TEMP*A(I,K)
  120                   CONTINUE
                     END IF
  140             CONTINUE
  160          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
C
            IF (UPPER) THEN
               DO 240 J = 1, N
                  DO 220 I = M, 1, -1
                     TEMP = B(I,J)
                     IF (NOCONJ) THEN
                        IF (NOUNIT) TEMP = TEMP*A(I,I)
                        DO 180 K = 1, I - 1
                           TEMP = TEMP + A(K,I)*B(K,J)
  180                   CONTINUE
                     ELSE
                        IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                        DO 200 K = 1, I - 1
                           TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  200                   CONTINUE
                     END IF
                     B(I,J) = ALPHA*TEMP
  220             CONTINUE
  240          CONTINUE
            ELSE
               DO 320 J = 1, N
                  DO 300 I = 1, M
                     TEMP = B(I,J)
                     IF (NOCONJ) THEN
                        IF (NOUNIT) TEMP = TEMP*A(I,I)
                        DO 260 K = I + 1, M
                           TEMP = TEMP + A(K,I)*B(K,J)
  260                   CONTINUE
                     ELSE
                        IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                        DO 280 K = I + 1, M
                           TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  280                   CONTINUE
                     END IF
                     B(I,J) = ALPHA*TEMP
  300             CONTINUE
  320          CONTINUE
            END IF
         END IF
      ELSE
         IF ((TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')) THEN
C
C           Form  B := alpha*B*A.
C
            IF (UPPER) THEN
               DO 400 J = N, 1, -1
                  TEMP = ALPHA
                  IF (NOUNIT) TEMP = TEMP*A(J,J)
                  DO 340 I = 1, M
                     B(I,J) = TEMP*B(I,J)
  340             CONTINUE
                  DO 380 K = 1, J - 1
                     IF (A(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*A(K,J)
                        DO 360 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  360                   CONTINUE
                     END IF
  380             CONTINUE
  400          CONTINUE
            ELSE
               DO 480 J = 1, N
                  TEMP = ALPHA
                  IF (NOUNIT) TEMP = TEMP*A(J,J)
                  DO 420 I = 1, M
                     B(I,J) = TEMP*B(I,J)
  420             CONTINUE
                  DO 460 K = J + 1, N
                     IF (A(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*A(K,J)
                        DO 440 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  440                   CONTINUE
                     END IF
  460             CONTINUE
  480          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
C
            IF (UPPER) THEN
               DO 560 K = 1, N
                  DO 520 J = 1, K - 1
                     IF (A(J,K).NE.ZERO) THEN
                        IF (NOCONJ) THEN
                           TEMP = ALPHA*A(J,K)
                        ELSE
                           TEMP = ALPHA*DCONJG(A(J,K))
                        END IF
                        DO 500 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  500                   CONTINUE
                     END IF
  520             CONTINUE
                  TEMP = ALPHA
                  IF (NOUNIT) THEN
                     IF (NOCONJ) THEN
                        TEMP = TEMP*A(K,K)
                     ELSE
                        TEMP = TEMP*DCONJG(A(K,K))
                     END IF
                  END IF
                  IF (TEMP.NE.ONE) THEN
                     DO 540 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  540                CONTINUE
                  END IF
  560          CONTINUE
            ELSE
               DO 640 K = N, 1, -1
                  DO 600 J = K + 1, N
                     IF (A(J,K).NE.ZERO) THEN
                        IF (NOCONJ) THEN
                           TEMP = ALPHA*A(J,K)
                        ELSE
                           TEMP = ALPHA*DCONJG(A(J,K))
                        END IF
                        DO 580 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  580                   CONTINUE
                     END IF
  600             CONTINUE
                  TEMP = ALPHA
                  IF (NOUNIT) THEN
                     IF (NOCONJ) THEN
                        TEMP = TEMP*A(K,K)
                     ELSE
                        TEMP = TEMP*DCONJG(A(K,K))
                     END IF
                  END IF
                  IF (TEMP.NE.ONE) THEN
                     DO 620 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  620                CONTINUE
                  END IF
  640          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06ZFF (ZTRMM ).
C
      END
