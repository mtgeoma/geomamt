      SUBROUTINE F06ZJF(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZTRSM  solves one of the matrix equations
C
C     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
C
C  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
C  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C
C     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
C
C  The matrix X is overwritten on B.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C           On entry, SIDE specifies whether op( A ) appears on the left
C           or right of X as follows:
C
C              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
C
C              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
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
C           contain  the  right-hand  side  matrix  B,  and  on exit  is
C           overwritten by the solution matrix  X.
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
      ENTRY             ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,
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
         CALL F06AAZ('F06ZJF/ZTRSM ',INFO)
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
C           Form  B := alpha*inv( A )*B.
C
            IF (UPPER) THEN
               DO 120 J = 1, N
                  IF (ALPHA.NE.ONE) THEN
                     DO 60 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
   60                CONTINUE
                  END IF
                  DO 100 K = M, 1, -1
                     IF (B(K,J).NE.ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        DO 80 I = 1, K - 1
                           B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                   CONTINUE
                     END IF
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 200 J = 1, N
                  IF (ALPHA.NE.ONE) THEN
                     DO 140 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
  140                CONTINUE
                  END IF
                  DO 180 K = 1, M
                     IF (B(K,J).NE.ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        DO 160 I = K + 1, M
                           B(I,J) = B(I,J) - B(K,J)*A(I,K)
  160                   CONTINUE
                     END IF
  180             CONTINUE
  200          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*inv( A' )*B
C           or    B := alpha*inv( conjg( A' ) )*B.
C
            IF (UPPER) THEN
               DO 280 J = 1, N
                  DO 260 I = 1, M
                     TEMP = ALPHA*B(I,J)
                     IF (NOCONJ) THEN
                        DO 220 K = 1, I - 1
                           TEMP = TEMP - A(K,I)*B(K,J)
  220                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(I,I)
                     ELSE
                        DO 240 K = 1, I - 1
                           TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  240                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                     END IF
                     B(I,J) = TEMP
  260             CONTINUE
  280          CONTINUE
            ELSE
               DO 360 J = 1, N
                  DO 340 I = M, 1, -1
                     TEMP = ALPHA*B(I,J)
                     IF (NOCONJ) THEN
                        DO 300 K = I + 1, M
                           TEMP = TEMP - A(K,I)*B(K,J)
  300                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(I,I)
                     ELSE
                        DO 320 K = I + 1, M
                           TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  320                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                     END IF
                     B(I,J) = TEMP
  340             CONTINUE
  360          CONTINUE
            END IF
         END IF
      ELSE
         IF ((TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')) THEN
C
C           Form  B := alpha*B*inv( A ).
C
            IF (UPPER) THEN
               DO 460 J = 1, N
                  IF (ALPHA.NE.ONE) THEN
                     DO 380 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
  380                CONTINUE
                  END IF
                  DO 420 K = 1, J - 1
                     IF (A(K,J).NE.ZERO) THEN
                        DO 400 I = 1, M
                           B(I,J) = B(I,J) - A(K,J)*B(I,K)
  400                   CONTINUE
                     END IF
  420             CONTINUE
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     DO 440 I = 1, M
                        B(I,J) = TEMP*B(I,J)
  440                CONTINUE
                  END IF
  460          CONTINUE
            ELSE
               DO 560 J = N, 1, -1
                  IF (ALPHA.NE.ONE) THEN
                     DO 480 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
  480                CONTINUE
                  END IF
                  DO 520 K = J + 1, N
                     IF (A(K,J).NE.ZERO) THEN
                        DO 500 I = 1, M
                           B(I,J) = B(I,J) - A(K,J)*B(I,K)
  500                   CONTINUE
                     END IF
  520             CONTINUE
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     DO 540 I = 1, M
                        B(I,J) = TEMP*B(I,J)
  540                CONTINUE
                  END IF
  560          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*inv( A' )
C           or    B := alpha*B*inv( conjg( A' ) ).
C
            IF (UPPER) THEN
               DO 660 K = N, 1, -1
                  IF (NOUNIT) THEN
                     IF (NOCONJ) THEN
                        TEMP = ONE/A(K,K)
                     ELSE
                        TEMP = ONE/DCONJG(A(K,K))
                     END IF
                     DO 580 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  580                CONTINUE
                  END IF
                  DO 620 J = 1, K - 1
                     IF (A(J,K).NE.ZERO) THEN
                        IF (NOCONJ) THEN
                           TEMP = A(J,K)
                        ELSE
                           TEMP = DCONJG(A(J,K))
                        END IF
                        DO 600 I = 1, M
                           B(I,J) = B(I,J) - TEMP*B(I,K)
  600                   CONTINUE
                     END IF
  620             CONTINUE
                  IF (ALPHA.NE.ONE) THEN
                     DO 640 I = 1, M
                        B(I,K) = ALPHA*B(I,K)
  640                CONTINUE
                  END IF
  660          CONTINUE
            ELSE
               DO 760 K = 1, N
                  IF (NOUNIT) THEN
                     IF (NOCONJ) THEN
                        TEMP = ONE/A(K,K)
                     ELSE
                        TEMP = ONE/DCONJG(A(K,K))
                     END IF
                     DO 680 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  680                CONTINUE
                  END IF
                  DO 720 J = K + 1, N
                     IF (A(J,K).NE.ZERO) THEN
                        IF (NOCONJ) THEN
                           TEMP = A(J,K)
                        ELSE
                           TEMP = DCONJG(A(J,K))
                        END IF
                        DO 700 I = 1, M
                           B(I,J) = B(I,J) - TEMP*B(I,K)
  700                   CONTINUE
                     END IF
  720             CONTINUE
                  IF (ALPHA.NE.ONE) THEN
                     DO 740 I = 1, M
                        B(I,K) = ALPHA*B(I,K)
  740                CONTINUE
                  END IF
  760          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06ZJF (ZTRSM ).
C
      END
