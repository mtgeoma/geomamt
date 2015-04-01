      SUBROUTINE F04LHX(SIDE,UPLO,TRANS,DIAG,M,N,M1,M2,A,LDA,B,LDB)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Purpose
C     =======
C
C     F04LHX  solves one of the matrix equations
C
C     op( A )*X = B,   or   X*op( A ) = B,
C
C     where X and B are m by n matrices, A is a unit, or non-unit, upper
C     or lower triangular matrix of the special forms
C
C     (  I  A12  A13 )  (upper)  or  (  I            )  (lower)
C     (     A22  A23 )               ( A21  A22      )
C     (           I  )               ( A31  A32   I  )
C
C     where the submatrices A22 in rows and columns m1 to m2 are
C     upper or lower triangular and the off-diagonal submatrices
C     are rectangular, and  op( A )  is one of
C
C     op( A ) = A   or   op( A ) = A'.
C
C     The matrix X is overwritten on B.
C
C     Parameters
C     ==========
C
C     SIDE   - CHARACTER*1.
C           On entry, SIDE specifies whether op( A ) appears on the left
C           or right of X as follows:
C
C              SIDE = 'L' or 'l'   op( A )*X = B.
C
C              SIDE = 'R' or 'r'   X*op( A ) = B.
C
C           Unchanged on exit.
C
C     UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix A is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C     TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANS = 'N' or 'n'   op( A ) = A.
C
C              TRANS = 'T' or 't'   op( A ) = A'.
C
C              TRANS = 'C' or 'c'   op( A ) = A'.
C
C           Unchanged on exit.
C
C     DIAG   - CHARACTER*1.
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
C     M      - INTEGER.
C           On entry, M specifies the number of rows of B. M must be at
C           least zero.
C           Unchanged on exit.
C
C     N      - INTEGER.
C           On entry, N specifies the number of columns of B.  N must be
C           at least zero.
C           Unchanged on exit.
C
C     M1     - INTEGER
C           On entry, M1 specifies the first row and column of the
C           diagonal submatrix A22. M1 must be at least 1.
C           Unchanged on exit.
C
C     M2     - INTEGER
C           On entry, M2 specifies the last row and column of the
C           diagonal submatrix A22. M2 must be at least m1-1, and must
C           not exceed k, where k is m when SIDE = 'L' or 'l' and is n
C           when SIDE = 'R'or 'r'.
C           Unchanged on exit.
C
C     A      - REAL        array of DIMENSION ( LDA, k ), where k is m
C           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C           Before entry  with  UPLO = 'U' or 'u', rows 1 to m2 and
C           columns m1 to k of the upper triangular part of the array
C           A must contain the elements of the submatrices A12, A13,
C           A22 and A23. The rest of the array is not referenced.
C           Before entry  with  UPLO = 'L' or 'l', rows m1 to k and
C           columns 1 to m2 of the lower triangular part of the array
C           A must contain the elements of the submatrices A21, A22,
C           A31 and A32. The rest of the array is not referenced.
C           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C           A22 are not referenced either,  but are assumed to be  unity
C           Unchanged on exit.
C
C     LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C           then LDA must be at least max( 1, n ).
C           Unchanged on exit.
C
C     B      - REAL             array of DIMENSION ( LDB, n ).
C           Before entry,  the leading  m by n part of the array  B must
C           contain  the  right-hand  side  matrix  B,  and  on exit  is
C           overwritten by the solution matrix  X.
C
C     LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   LDB  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C     Level 3 Blas routine.
C
C     -- Written on 20-February-1987.
C     Sven Hammarling, Nag Central Office.
C     Richard W. Brankin, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER           LDA, LDB, M, M1, M2, N
      CHARACTER*1       DIAG, SIDE, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*)
C     .. Local Scalars ..
      INTEGER           INFO, J, NROWA
      LOGICAL           LSIDE
      CHARACTER*1       NOTR
C     .. External Subroutines ..
      EXTERNAL          F04LHW, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      LSIDE = SIDE .EQ. 'L' .OR. SIDE .EQ. 'l'
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      INFO = 0
      IF (( .NOT. LSIDE) .AND. (SIDE.NE.'R' .AND. SIDE.NE.'r')) THEN
         INFO = 1
      ELSE IF ((UPLO.NE.'U' .AND. UPLO.NE.'u')
     *         .AND. (UPLO.NE.'L' .AND. UPLO.NE.'l')) THEN
         INFO = 2
      ELSE IF ((TRANS.NE.'N' .AND. TRANS.NE.'n')
     *         .AND. (TRANS.NE.'T' .AND. TRANS.NE.'t')
     *         .AND. (TRANS.NE.'C' .AND. TRANS.NE.'c')) THEN
         INFO = 3
      ELSE IF ((DIAG.NE.'U' .AND. DIAG.NE.'u')
     *         .AND. (DIAG.NE.'N' .AND. DIAG.NE.'n')) THEN
         INFO = 4
      ELSE IF (M.LT.0) THEN
         INFO = 5
      ELSE IF (N.LT.0) THEN
         INFO = 6
      ELSE IF (M1.LT.1) THEN
         INFO = 7
      ELSE IF ((M2.LT.(M1-1)) .OR. (LSIDE .AND. M2.GT.M)
     *         .OR. ( .NOT. LSIDE .AND. M2.GT.N)) THEN
         INFO = 8
      ELSE IF (((UPLO.EQ.'U' .OR. UPLO.EQ.'u') .AND. LDA.LT.M2)
     *         .OR. ((UPLO.EQ.'L' .OR. UPLO.EQ.'l') .AND. LDA.LT.MAX(1,
     *         NROWA))) THEN
         INFO = 10
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F04LHX       ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF (N.EQ.0) RETURN
C
C     Start the operations.
C
      IF (LSIDE) THEN
C
C        Solve  op( A )*X = B.
C
         DO 20 J = 1, N
            CALL F04LHW(UPLO,TRANS,DIAG,M,M1,M2,A,LDA,B(1,J),1)
   20    CONTINUE
      ELSE
         IF (TRANS.EQ.'N' .OR. TRANS.EQ.'n') THEN
            NOTR = 'T'
         ELSE
            NOTR = 'N'
         END IF
C
C        Solve  X*op( A ) = B  using  op( A' )*X' = B'.
C
         DO 40 J = 1, M
            CALL F04LHW(UPLO,NOTR,DIAG,N,M1,M2,A,LDA,B(J,1),LDB)
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F04LHX .
C
      END
