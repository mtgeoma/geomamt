      SUBROUTINE F06ZRF(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1647 (JUN 1995).
C
C  Purpose
C  =======
C
C  ZHER2K  performs one of the hermitian rank 2k operations
C
C     C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) + beta*C,
C
C  or
C
C     C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A + beta*C,
C
C  where  alpha and beta  are scalars with  beta  real,  C is an  n by n
C  hermitian matrix and  A and B  are  n by k matrices in the first case
C  and  k by n  matrices in the second case.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On  entry,   UPLO  specifies  whether  the  upper  or  lower
C           triangular  part  of the  array  C  is to be  referenced  as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry,  TRANS  specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'    C := alpha*A*conjg( B' )          +
C                                         conjg( alpha )*B*conjg( A' ) +
C                                         beta*C.
C
C              TRANS = 'C' or 'c'    C := alpha*conjg( A' )*B          +
C                                         conjg( alpha )*conjg( B' )*A +
C                                         beta*C.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N specifies the order of the matrix C.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
C           of  columns  of the  matrices  A and B,  and on  entry  with
C           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
C           matrices  A and B.  K must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by n  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
C           then  LDA must be at least  max( 1, n ), otherwise  LDA must
C           be at least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  k by n  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
C           then  LDB must be at least  max( 1, n ), otherwise  LDB must
C           be at least  max( 1, k ).
C           Unchanged on exit.
C
C  BETA   - REAL            .
C           On entry, BETA specifies the scalar beta.
C           Unchanged on exit.
C
C  C      - COMPLEX          array of DIMENSION ( LDC, n ).
C           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
C           upper triangular part of the array C must contain the upper
C           triangular part  of the  hermitian matrix  and the strictly
C           lower triangular part of C is not referenced.  On exit, the
C           upper triangular part of the array  C is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
C           lower triangular part of the array C must contain the lower
C           triangular part  of the  hermitian matrix  and the strictly
C           upper triangular part of C is not referenced.  On exit, the
C           lower triangular part of the array  C is overwritten by the
C           lower triangular part of the updated matrix.
C           Note that the imaginary parts of the diagonal elements need
C           not be set,  they are assumed to be zero,  and on exit they
C           are set to zero.
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, n ).
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
C  -- Modified on 28-March-1995 to set C(J,J) to REAL( C(J,J) ) when
C     BETA = 1.
C     Sven Hammarling (following the correction by Ed Anderson).
C
C     .. Entry Points ..
      ENTRY             ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA
      DOUBLE PRECISION  BETA
      INTEGER           K, LDA, LDB, LDC, N
      CHARACTER*1       TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP1, TEMP2
      INTEGER           I, INFO, J, L, NROWA
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCONJG, MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      IF ((TRANS.EQ.'N' .OR. TRANS.EQ.'n')) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
      INFO = 0
      IF (( .NOT. UPPER) .AND. ( .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')))
     *    THEN
         INFO = 1
      ELSE IF (( .NOT. (TRANS.EQ.'N' .OR. TRANS.EQ.'n'))
     *         .AND. ( .NOT. (TRANS.EQ.'C' .OR. TRANS.EQ.'c'))) THEN
         INFO = 2
      ELSE IF (N.LT.0) THEN
         INFO = 3
      ELSE IF (K.LT.0) THEN
         INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
         INFO = 12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZRF/ZHER2K',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *    .AND. (BETA.EQ.ONE))) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         IF (UPPER) THEN
            IF (BETA.EQ.DBLE(ZERO)) THEN
               DO 40 J = 1, N
                  DO 20 I = 1, J
                     C(I,J) = ZERO
   20             CONTINUE
   40          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 60 I = 1, J - 1
                     C(I,J) = BETA*C(I,J)
   60             CONTINUE
                  C(J,J) = BETA*DBLE(C(J,J))
   80          CONTINUE
            END IF
         ELSE
            IF (BETA.EQ.DBLE(ZERO)) THEN
               DO 120 J = 1, N
                  DO 100 I = J, N
                     C(I,J) = ZERO
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 160 J = 1, N
                  C(J,J) = BETA*DBLE(C(J,J))
                  DO 140 I = J + 1, N
                     C(I,J) = BETA*C(I,J)
  140             CONTINUE
  160          CONTINUE
            END IF
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF ((TRANS.EQ.'N' .OR. TRANS.EQ.'n')) THEN
C
C        Form  C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) +
C                   C.
C
         IF (UPPER) THEN
            DO 260 J = 1, N
               IF (BETA.EQ.DBLE(ZERO)) THEN
                  DO 180 I = 1, J
                     C(I,J) = ZERO
  180             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 200 I = 1, J - 1
                     C(I,J) = BETA*C(I,J)
  200             CONTINUE
                  C(J,J) = BETA*DBLE(C(J,J))
               ELSE
                  C(J,J) = DBLE(C(J,J))
               END IF
               DO 240 L = 1, K
                  IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                     TEMP1 = ALPHA*DCONJG(B(J,L))
                     TEMP2 = DCONJG(ALPHA*A(J,L))
                     DO 220 I = 1, J - 1
                        C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  220                CONTINUE
                     C(J,J) = DBLE(C(J,J)) + DBLE(A(J,L)*TEMP1+B(J,L)
     *                        *TEMP2)
                  END IF
  240          CONTINUE
  260       CONTINUE
         ELSE
            DO 360 J = 1, N
               IF (BETA.EQ.DBLE(ZERO)) THEN
                  DO 280 I = J, N
                     C(I,J) = ZERO
  280             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 300 I = J + 1, N
                     C(I,J) = BETA*C(I,J)
  300             CONTINUE
                  C(J,J) = BETA*DBLE(C(J,J))
               ELSE
                  C(J,J) = DBLE(C(J,J))
               END IF
               DO 340 L = 1, K
                  IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                     TEMP1 = ALPHA*DCONJG(B(J,L))
                     TEMP2 = DCONJG(ALPHA*A(J,L))
                     DO 320 I = J + 1, N
                        C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  320                CONTINUE
                     C(J,J) = DBLE(C(J,J)) + DBLE(A(J,L)*TEMP1+B(J,L)
     *                        *TEMP2)
                  END IF
  340          CONTINUE
  360       CONTINUE
         END IF
      ELSE
C
C        Form  C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A +
C                   C.
C
         IF (UPPER) THEN
            DO 420 J = 1, N
               DO 400 I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 380 L = 1, K
                     TEMP1 = TEMP1 + DCONJG(A(L,I))*B(L,J)
                     TEMP2 = TEMP2 + DCONJG(B(L,I))*A(L,J)
  380             CONTINUE
                  IF (I.EQ.J) THEN
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(J,J) = DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     ELSE
                        C(J,J) = BETA*DBLE(C(J,J)) +
     *                           DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     END IF
                  ELSE
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(I,J) = ALPHA*TEMP1 + DCONJG(ALPHA)*TEMP2
                     ELSE
                        C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     *                           DCONJG(ALPHA)*TEMP2
                     END IF
                  END IF
  400          CONTINUE
  420       CONTINUE
         ELSE
            DO 480 J = 1, N
               DO 460 I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 440 L = 1, K
                     TEMP1 = TEMP1 + DCONJG(A(L,I))*B(L,J)
                     TEMP2 = TEMP2 + DCONJG(B(L,I))*A(L,J)
  440             CONTINUE
                  IF (I.EQ.J) THEN
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(J,J) = DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     ELSE
                        C(J,J) = BETA*DBLE(C(J,J)) +
     *                           DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     END IF
                  ELSE
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(I,J) = ALPHA*TEMP1 + DCONJG(ALPHA)*TEMP2
                     ELSE
                        C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     *                           DCONJG(ALPHA)*TEMP2
                     END IF
                  END IF
  460          CONTINUE
  480       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06ZRF (ZHER2K).
C
      END
