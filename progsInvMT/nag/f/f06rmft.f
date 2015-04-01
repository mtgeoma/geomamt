      DOUBLE PRECISION FUNCTION F06RMF(NORM,N,A,LDA,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 DLANHS
C     ENTRY                DLANHS(NORM,N,A,LDA,WORK)
C
C  Purpose
C  =======
C
C  DLANHS  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the  element of  largest absolute value  of a
C  Hessenberg matrix A.
C
C  Description
C  ===========
C
C  DLANHS returns the value
C
C     DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
C          Specifies the value to be returned in DLANHS as described
C          above.
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
C          set to zero.
C
C  A       (input) REAL array, dimension (LDA,N)
C          The n by n upper Hessenberg matrix A; the part of A below the
C          first sub-diagonal is not referenced.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(N,1).
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
      INTEGER                          LDA, N
      CHARACTER                        NORM
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(LDA,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUM, VALUE
      INTEGER                          I, J
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
            DO 20 I = 1, MIN(N,J+1)
               VALUE = MAX(VALUE,ABS(A(I,J)))
   20       CONTINUE
   40    CONTINUE
      ELSE IF (((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1')) THEN
C
C        Find norm1(A).
C
         VALUE = ZERO
         DO 80 J = 1, N
            SUM = ZERO
            DO 60 I = 1, MIN(N,J+1)
               SUM = SUM + ABS(A(I,J))
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
            DO 120 I = 1, MIN(N,J+1)
               WORK(I) = WORK(I) + ABS(A(I,J))
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
            CALL F06FJF(MIN(N,J+1),A(1,J),1,SCALE,SUM)
  180    CONTINUE
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06RMF = VALUE
      RETURN
C
C     End of F06RMF (DLANHS)
C
      END
