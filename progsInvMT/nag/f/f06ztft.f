      SUBROUTINE F06ZTF(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZSYMM  performs one of the matrix-matrix operations
C
C     C := alpha*A*B + beta*C,
C
C  or
C
C     C := alpha*B*A + beta*C,
C
C  where  alpha and beta are scalars, A is a symmetric matrix and  B and
C  C are m by n matrices.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C           On entry,  SIDE  specifies whether  the  symmetric matrix  A
C           appears on the  left or right  in the  operation as follows:
C
C              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
C
C              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
C
C           Unchanged on exit.
C
C  UPLO   - CHARACTER*1.
C           On  entry,   UPLO  specifies  whether  the  upper  or  lower
C           triangular  part  of  the  symmetric  matrix   A  is  to  be
C           referenced as follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of the
C                                  symmetric matrix is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of the
C                                  symmetric matrix is to be referenced.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies the number of rows of the matrix  C.
C           M  must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix C.
C           N  must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
C           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
C           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
C           the array  A  must contain the  symmetric matrix,  such that
C           when  UPLO = 'U' or 'u', the leading m by m upper triangular
C           part of the array  A  must contain the upper triangular part
C           of the  symmetric matrix and the  strictly  lower triangular
C           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
C           the leading  m by m  lower triangular part  of the  array  A
C           must  contain  the  lower triangular part  of the  symmetric
C           matrix and the  strictly upper triangular part of  A  is not
C           referenced.
C           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
C           the array  A  must contain the  symmetric matrix,  such that
C           when  UPLO = 'U' or 'u', the leading n by n upper triangular
C           part of the array  A  must contain the upper triangular part
C           of the  symmetric matrix and the  strictly  lower triangular
C           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
C           the leading  n by n  lower triangular part  of the  array  A
C           must  contain  the  lower triangular part  of the  symmetric
C           matrix and the  strictly upper triangular part of  A  is not
C           referenced.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least max( 1, n ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, n ).
C           Before entry, the leading  m by n part of the array  B  must
C           contain the matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   LDB  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C  BETA   - COMPLEX         .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - COMPLEX          array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n updated
C           matrix.
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
      ENTRY             ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           LDA, LDB, LDC, M, N
      CHARACTER*1       SIDE, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP1, TEMP2
      INTEGER           I, INFO, J, K, NROWA
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Set NROWA as the number of rows of A.
C
      IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
C     Test the input parameters.
C
      INFO = 0
      IF (( .NOT. (SIDE.EQ.'L' .OR. SIDE.EQ.'l'))
     *    .AND. ( .NOT. (SIDE.EQ.'R' .OR. SIDE.EQ.'r'))) THEN
         INFO = 1
      ELSE IF (( .NOT. UPPER) .AND. ( .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.
     *         'l'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 7
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 9
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZTF/ZSYMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. ((ALPHA.EQ.ZERO)
     *    .AND. (BETA.EQ.ONE))) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  C(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
         ELSE
            DO 80 J = 1, N
               DO 60 I = 1, M
                  C(I,J) = BETA*C(I,J)
   60          CONTINUE
   80       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C        Form  C := alpha*A*B + beta*C.
C
         IF (UPPER) THEN
            DO 140 J = 1, N
               DO 120 I = 1, M
                  TEMP1 = ALPHA*B(I,J)
                  TEMP2 = ZERO
                  DO 100 K = 1, I - 1
                     C(K,J) = C(K,J) + TEMP1*A(K,I)
                     TEMP2 = TEMP2 + B(K,J)*A(K,I)
  100             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                  ELSE
                     C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2
                  END IF
  120          CONTINUE
  140       CONTINUE
         ELSE
            DO 200 J = 1, N
               DO 180 I = M, 1, -1
                  TEMP1 = ALPHA*B(I,J)
                  TEMP2 = ZERO
                  DO 160 K = I + 1, M
                     C(K,J) = C(K,J) + TEMP1*A(K,I)
                     TEMP2 = TEMP2 + B(K,J)*A(K,I)
  160             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                  ELSE
                     C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2
                  END IF
  180          CONTINUE
  200       CONTINUE
         END IF
      ELSE
C
C        Form  C := alpha*B*A + beta*C.
C
         DO 340 J = 1, N
            TEMP1 = ALPHA*A(J,J)
            IF (BETA.EQ.ZERO) THEN
               DO 220 I = 1, M
                  C(I,J) = TEMP1*B(I,J)
  220          CONTINUE
            ELSE
               DO 240 I = 1, M
                  C(I,J) = BETA*C(I,J) + TEMP1*B(I,J)
  240          CONTINUE
            END IF
            DO 280 K = 1, J - 1
               IF (UPPER) THEN
                  TEMP1 = ALPHA*A(K,J)
               ELSE
                  TEMP1 = ALPHA*A(J,K)
               END IF
               DO 260 I = 1, M
                  C(I,J) = C(I,J) + TEMP1*B(I,K)
  260          CONTINUE
  280       CONTINUE
            DO 320 K = J + 1, N
               IF (UPPER) THEN
                  TEMP1 = ALPHA*A(J,K)
               ELSE
                  TEMP1 = ALPHA*A(K,J)
               END IF
               DO 300 I = 1, M
                  C(I,J) = C(I,J) + TEMP1*B(I,K)
  300          CONTINUE
  320       CONTINUE
  340    CONTINUE
      END IF
C
      RETURN
C
C     End of F06ZTF (ZSYMM ).
C
      END
