      SUBROUTINE F06YPF(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  DSYRK  performs one of the symmetric rank k operations
C
C     C := alpha*A*A' + beta*C,
C
C  or
C
C     C := alpha*A'*A + beta*C,
C
C  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
C  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
C  in the second case.
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
C              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
C
C              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
C
C              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
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
C           of  columns   of  the   matrix   A,   and  on   entry   with
C           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
C           of rows of the matrix  A.  K must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - REAL            .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
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
C  BETA   - REAL            .
C           On entry, BETA specifies the scalar beta.
C           Unchanged on exit.
C
C  C      - REAL             array of DIMENSION ( LDC, n ).
C           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
C           upper triangular part of the array C must contain the upper
C           triangular part  of the  symmetric matrix  and the strictly
C           lower triangular part of C is not referenced.  On exit, the
C           upper triangular part of the array  C is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
C           lower triangular part of the array C must contain the lower
C           triangular part  of the  symmetric matrix  and the strictly
C           upper triangular part of C is not referenced.  On exit, the
C           lower triangular part of the array  C is overwritten by the
C           lower triangular part of the updated matrix.
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
C
C     .. Entry Points ..
      ENTRY             DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           K, LDA, LDC, N
      CHARACTER*1       TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(LDC,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, INFO, J, L, NROWA
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
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
     *         .AND. ( .NOT. (TRANS.EQ.'T' .OR. TRANS.EQ.'t'))
     *         .AND. ( .NOT. (TRANS.EQ.'C' .OR. TRANS.EQ.'c'))) THEN
         INFO = 2
      ELSE IF (N.LT.0) THEN
         INFO = 3
      ELSE IF (K.LT.0) THEN
         INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
         INFO = 10
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06YPF/DSYRK ',INFO)
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
            IF (BETA.EQ.ZERO) THEN
               DO 40 J = 1, N
                  DO 20 I = 1, J
                     C(I,J) = ZERO
   20             CONTINUE
   40          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 60 I = 1, J
                     C(I,J) = BETA*C(I,J)
   60             CONTINUE
   80          CONTINUE
            END IF
         ELSE
            IF (BETA.EQ.ZERO) THEN
               DO 120 J = 1, N
                  DO 100 I = J, N
                     C(I,J) = ZERO
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 140 I = J, N
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
C        Form  C := alpha*A*A' + beta*C.
C
         IF (UPPER) THEN
            DO 260 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 180 I = 1, J
                     C(I,J) = ZERO
  180             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 200 I = 1, J
                     C(I,J) = BETA*C(I,J)
  200             CONTINUE
               END IF
               DO 240 L = 1, K
                  IF (A(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*A(J,L)
                     DO 220 I = 1, J
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  220                CONTINUE
                  END IF
  240          CONTINUE
  260       CONTINUE
         ELSE
            DO 360 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 280 I = J, N
                     C(I,J) = ZERO
  280             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 300 I = J, N
                     C(I,J) = BETA*C(I,J)
  300             CONTINUE
               END IF
               DO 340 L = 1, K
                  IF (A(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*A(J,L)
                     DO 320 I = J, N
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  320                CONTINUE
                  END IF
  340          CONTINUE
  360       CONTINUE
         END IF
      ELSE
C
C        Form  C := alpha*A'*A + beta*C.
C
         IF (UPPER) THEN
            DO 420 J = 1, N
               DO 400 I = 1, J
                  TEMP = ZERO
                  DO 380 L = 1, K
                     TEMP = TEMP + A(L,I)*A(L,J)
  380             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  400          CONTINUE
  420       CONTINUE
         ELSE
            DO 480 J = 1, N
               DO 460 I = J, N
                  TEMP = ZERO
                  DO 440 L = 1, K
                     TEMP = TEMP + A(L,I)*A(L,J)
  440             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  460          CONTINUE
  480       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06YPF (DSYRK ).
C
      END
