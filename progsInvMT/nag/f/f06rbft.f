      DOUBLE PRECISION FUNCTION F06RBF(NORM,N,KL,KU,AB,LDAB,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 DLANGB
C     ENTRY                DLANGB(NORM,N,KL,KU,AB,LDAB,WORK)
C
C  Purpose
C  =======
C
C  DLANGB  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the element of  largest absolute value  of an
C  n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.
C
C  Description
C  ===========
C
C  DLANGB returns the value
C
C     DLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
C              (
C              ( norm1(A),         NORM = '1', 'O' or 'o'
C              (
C              ( normI(A),         NORM = 'I' or 'i'
C              (
C              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
C
C  where  norm1  denotes the  one norm of a matrix (maximum column sum),
C  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
C  normF  denotes the  Frobenius norm of a matrix (square root of sum of
C  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
C
C  Arguments
C  =========
C
C  NORM    (input) CHARACTER*1
C          Specifies the value to be returned in DLANGB as described
C          above.
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.  When N = 0, DLANGB is
C          set to zero.
C
C  KL      (input) INTEGER
C          The number of sub-diagonals of the matrix A.  KL >= 0.
C
C  KU      (input) INTEGER
C          The number of super-diagonals of the matrix A.  KU >= 0.
C
C  AB      (input) REAL array, dimension (LDAB,N)
C          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th
C          column of A is stored in the j-th column of the array AB as
C          follows:
C          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KL+KU+1.
C
C  WORK    (workspace) REAL array, dimension (LWORK),
C          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
C          referenced.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, ZERO
      PARAMETER                        (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER                          KL, KU, LDAB, N
      CHARACTER                        NORM
C     .. Array Arguments ..
      DOUBLE PRECISION                 AB(LDAB,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUM, VALUE
      INTEGER                          I, J, K, L
C     .. External Subroutines ..
      EXTERNAL                         F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      IF (N.EQ.0) THEN
         VALUE = ZERO
      ELSE IF ((NORM.EQ.'M' .OR. NORM.EQ.'m')) THEN
C
C        Find max(abs(A(i,j))).
C
         VALUE = ZERO
         DO 40 J = 1, N
            DO 20 I = MAX(KU+2-J,1), MIN(N+KU+1-J,KL+KU+1)
               VALUE = MAX(VALUE,ABS(AB(I,J)))
   20       CONTINUE
   40    CONTINUE
      ELSE IF (((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1')) THEN
C
C        Find norm1(A).
C
         VALUE = ZERO
         DO 80 J = 1, N
            SUM = ZERO
            DO 60 I = MAX(KU+2-J,1), MIN(N+KU+1-J,KL+KU+1)
               SUM = SUM + ABS(AB(I,J))
   60       CONTINUE
            VALUE = MAX(VALUE,SUM)
   80    CONTINUE
      ELSE IF ((NORM.EQ.'I' .OR. NORM.EQ.'i')) THEN
C
C        Find normI(A).
C
         DO 100 I = 1, N
            WORK(I) = ZERO
  100    CONTINUE
         DO 140 J = 1, N
            K = KU + 1 - J
            DO 120 I = MAX(1,J-KU), MIN(N,J+KL)
               WORK(I) = WORK(I) + ABS(AB(K+I,J))
  120       CONTINUE
  140    CONTINUE
         VALUE = ZERO
         DO 160 I = 1, N
            VALUE = MAX(VALUE,WORK(I))
  160    CONTINUE
      ELSE IF (((NORM.EQ.'F' .OR. NORM.EQ.'f'))
     *         .OR. ((NORM.EQ.'E' .OR. NORM.EQ.'e'))) THEN
C
C        Find normF(A).
C
         SCALE = ZERO
         SUM = ONE
         DO 180 J = 1, N
            L = MAX(1,J-KU)
            K = KU + 1 - J + L
            CALL F06FJF(MIN(N,J+KL)-L+1,AB(K,J),1,SCALE,SUM)
  180    CONTINUE
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06RBF = VALUE
      RETURN
C
C     End of F06RBF (DLANGB)
C
      END
