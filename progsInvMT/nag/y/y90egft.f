      DOUBLE PRECISION FUNCTION Y90EGF(NORM,MATRIX,M,N,A,LDA)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C  This routine originally called F06VGF.
C
C  Purpose
C  =======
C
C  Y90EGF  returns the value  of the  one norm,  or the  Frobenius norm,
C  or the  element of  largest absolute value  of a  complex  matrix  A.
C  A  may be  rectangular,  or  square,  or  triangular,  or  hermitian.
C
C  Description
C  ===========
C
C  Y90EGF returns the value
C
C     Y90EGF = ( max( abs( a( i, j ) ) ) , NORM = 'M' or 'm'
C              (
C              ( norm1( A )  ,             NORM = '1', 'O' or 'o'
C              (
C              ( normF( A ) ,              NORM = 'F', 'f', 'E' or 'e'
C
C  where norm1 denotes the one norm of a matrix (maximum column sum) and
C  normF  denotes the  Frobenius norm of a matrix (square root of sum of
C  squares).  Note that  max( abs( a( i, j ) ) )  is not a  matrix norm.
C
C  The  type of matrix  for which  Y90EGF  is returned  is determined by
C  the parameter MATRIX as follows.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded  as a general matrix,
C  If   MATRIX = 'U' or 'u'   then  A  is regarded  as upper triangular,
C  If   MATRIX = 'L' or 'l'   then  A  is regarded  as lower triangular,
C  If   MATRIX = 'H' or 'h'   then  A  is regarded as hermitian and only
C                             the  upper triangular part of the array  A
C                             is referenced,
C  If   MATRIX = 'E'  or 'e'  then  A  is regarded as hermitian and only
C                             the  lower triangular part of the array  A
C                             is referenced.
C
C  Parameters
C  ==========
C
C  NORM  -  CHARACTER*1.
C
C           On entry,  NORM specifies the value to be returned in Y90EGF
C           as described above.
C
C           Unchanged on exit.
C
C  MATRIX - CHARACTER*1.
C
C           On entry,  MATRIX  specifies the type of  matrix and, in the
C           case of an  hermitian matrix, the part of the array in which
C           the matrix is stored as described above.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies the number of rows of the matrix  A.
C           M  must be at least  zero and when the  matrix is  hermitian
C           then  M must be equal to N.  When  M = 0  then Y90EGF is set
C           to zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.  When  N = 0  then Y90EGF is set to
C           zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  A      - COMPLEX array of DIMENSION ( LDA, n ).
C
C           Before entry,  A  must contain the  m by n  matrix for which
C           Y90EGF is required.
C
C           If   MATRIX = 'U' or 'u' or 'H' or 'h'   then  the  strictly
C           lower triangular part of A is not referenced.
C
C           If   MATRIX = 'L' or 'l' or 'E' or 'e'   then  the  strictly
C           upper triangular part of A is not referenced.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On entry, LDA specifies the first dimension of A as declared
C           in  the  calling  (sub)  program.  LDA  must be at least  M.
C
C           Unchanged on exit.
C
C  Further comments
C  ================
C
C  If A is part of a matrix B partitioned as
C
C     B = ( B1  B2 ) ,
C         ( B3  A  )
C
C  where  B1 is an  l by k matrix ( l.ge.0, k.ge.0 ),  then this routine
C  may be called with the parameter  A as  b( l + 1, k + 1 ) and  LDA as
C  the first dimension of  B  as declared in the  calling (sub) program.
C
C  This routine can be  inefficient on paged machines  when the one norm
C  is required, the matrix is hermitian and N is large.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 14-November-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, ZERO
      PARAMETER                        (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER                          LDA, M, N
      CHARACTER*1                      MATRIX, NORM
C     .. Array Arguments ..
      COMPLEX*16                       A(LDA,*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUM, VALUE
      INTEGER                          I, J
C     .. External Functions ..
      DOUBLE PRECISION                 F06BMF
      EXTERNAL                         F06BMF
C     .. External Subroutines ..
      EXTERNAL                         F06KJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, MIN
C     .. Executable Statements ..
      IF (MIN(M,N).EQ.0) THEN
         VALUE = ZERO
      ELSE IF ((NORM.EQ.'M') .OR. (NORM.EQ.'m')) THEN
C
C        Find  max( abs( a( i, j ) ) ).
C
         VALUE = ZERO
         IF ((MATRIX.EQ.'G') .OR. (MATRIX.EQ.'g')) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  VALUE = MAX(VALUE,ABS(A(I,J)))
   20          CONTINUE
   40       CONTINUE
         ELSE IF ((MATRIX.EQ.'U') .OR. (MATRIX.EQ.'u')
     *            .OR. (MATRIX.EQ.'H') .OR. (MATRIX.EQ.'h')) THEN
            DO 80 J = 1, N
               DO 60 I = 1, MIN(M,J)
                  VALUE = MAX(VALUE,ABS(A(I,J)))
   60          CONTINUE
   80       CONTINUE
         ELSE IF ((MATRIX.EQ.'L') .OR. (MATRIX.EQ.'l')
     *            .OR. (MATRIX.EQ.'E') .OR. (MATRIX.EQ.'e')) THEN
            DO 120 J = 1, MIN(M,N)
               DO 100 I = J, M
                  VALUE = MAX(VALUE,ABS(A(I,J)))
  100          CONTINUE
  120       CONTINUE
         END IF
      ELSE IF ((NORM.EQ.'1') .OR. (NORM.EQ.'O') .OR. (NORM.EQ.'o')) THEN
C
C        Find  norm1( A ).
C
         VALUE = ZERO
         IF ((MATRIX.EQ.'G') .OR. (MATRIX.EQ.'g')) THEN
            DO 160 J = 1, N
               SUM = ZERO
               DO 140 I = 1, M
                  SUM = SUM + ABS(A(I,J))
  140          CONTINUE
               VALUE = MAX(VALUE,SUM)
  160       CONTINUE
         ELSE IF ((MATRIX.EQ.'U') .OR. (MATRIX.EQ.'u')) THEN
            DO 200 J = 1, N
               SUM = ZERO
               DO 180 I = 1, MIN(M,J)
                  SUM = SUM + ABS(A(I,J))
  180          CONTINUE
               VALUE = MAX(VALUE,SUM)
  200       CONTINUE
         ELSE IF ((MATRIX.EQ.'L') .OR. (MATRIX.EQ.'l')) THEN
            DO 240 J = 1, MIN(M,N)
               SUM = ZERO
               DO 220 I = J, M
                  SUM = SUM + ABS(A(I,J))
  220          CONTINUE
               VALUE = MAX(VALUE,SUM)
  240       CONTINUE
         ELSE IF ((MATRIX.EQ.'H') .OR. (MATRIX.EQ.'h')) THEN
            DO 300 J = 1, N
               SUM = ZERO
               DO 260 I = 1, J
                  SUM = SUM + ABS(A(I,J))
  260          CONTINUE
               DO 280 I = J + 1, N
                  SUM = SUM + ABS(A(J,I))
  280          CONTINUE
               VALUE = MAX(VALUE,SUM)
  300       CONTINUE
         ELSE IF ((MATRIX.EQ.'E') .OR. (MATRIX.EQ.'e')) THEN
            DO 360 J = 1, N
               SUM = ZERO
               DO 320 I = 1, J - 1
                  SUM = SUM + ABS(A(J,I))
  320          CONTINUE
               DO 340 I = J, N
                  SUM = SUM + ABS(A(I,J))
  340          CONTINUE
               VALUE = MAX(VALUE,SUM)
  360       CONTINUE
         END IF
      ELSE IF ((NORM.EQ.'F') .OR. (NORM.EQ.'f') .OR. (NORM.EQ.'E')
     *         .OR. (NORM.EQ.'e')) THEN
C
C        Find  normF( A ).
C
         SCALE = ZERO
         SUM = ONE
         IF ((MATRIX.EQ.'G') .OR. (MATRIX.EQ.'g')) THEN
            DO 380 J = 1, N
               CALL F06KJF(M,A(1,J),1,SCALE,SUM)
  380       CONTINUE
         ELSE IF ((MATRIX.EQ.'U') .OR. (MATRIX.EQ.'u')) THEN
            DO 400 J = 1, N
               CALL F06KJF(MIN(M,J),A(1,J),1,SCALE,SUM)
  400       CONTINUE
         ELSE IF ((MATRIX.EQ.'L') .OR. (MATRIX.EQ.'l')) THEN
            DO 420 J = 1, MIN(M,N)
               CALL F06KJF(M-J+1,A(J,J),1,SCALE,SUM)
  420       CONTINUE
         ELSE IF ((MATRIX.EQ.'H') .OR. (MATRIX.EQ.'h')
     *            .OR. (MATRIX.EQ.'E') .OR. (MATRIX.EQ.'e')) THEN
            IF ((MATRIX.EQ.'H') .OR. (MATRIX.EQ.'h')) THEN
               DO 440 J = 2, N
                  CALL F06KJF(J-1,A(1,J),1,SCALE,SUM)
  440          CONTINUE
            ELSE
               DO 460 J = 1, N - 1
                  CALL F06KJF(N-J,A(J+1,J),1,SCALE,SUM)
  460          CONTINUE
            END IF
            SUM = 2*SUM
            CALL F06KJF(N,A(1,1),LDA+1,SCALE,SUM)
         END IF
         VALUE = F06BMF(SCALE,SUM)
      END IF
C
      Y90EGF = VALUE
      RETURN
C
C     End of Y90EGF.
C
      END
