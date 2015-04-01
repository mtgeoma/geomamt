      SUBROUTINE F07HSF(UPLO,N,KD,NRHS,AB,LDAB,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZPBTRS(UPLO,N,KD,NRHS,AB,LDAB,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  ZPBTRS solves a system of linear equations A*X = B with a Hermitian
C  positive definite band matrix A using the Cholesky factorization
C  A = U'*U or A = L*L' computed by F07HRF.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the factor stored in A is upper or lower
C          triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  KD      (input) INTEGER
C          The number of super-diagonals of the matrix A if UPLO = 'U',
C          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
C
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrix B.  NRHS >= 0.
C
C  AB      (input) COMPLEX array, dimension (LDAB,N)
C          The triangular factor U or L from the Cholesky factorization
C          A = U'*U or A = L*L' of the band matrix A, stored in the
C          first KD+1 rows of the array.  The j-th column of U or L is
C          stored in the array AB as follows:
C          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
C          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KD+1.
C
C  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
C          On entry, the right hand side vectors B for the system of
C          linear equations.
C          On exit, the solution vectors, X.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, KD, LDAB, LDB, N, NRHS
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AB(LDAB,*), B(LDB,*)
C     .. Local Scalars ..
      INTEGER           J
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          ZTBSV, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (KD.LT.0) THEN
         INFO = -3
      ELSE IF (NRHS.LT.0) THEN
         INFO = -4
      ELSE IF (LDAB.LT.KD+1) THEN
         INFO = -6
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07HSF/ZPBTRS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Solve A*X = B where A = U'*U.
C
         DO 20 J = 1, NRHS
C
C           Solve U'*X = B, overwriting B with X.
C
            CALL ZTBSV('Upper','Conjugate transpose','Non-unit',N,KD,AB,
     *                 LDAB,B(1,J),1)
C
C           Solve U*X = B, overwriting B with X.
C
            CALL ZTBSV('Upper','No transpose','Non-unit',N,KD,AB,LDAB,
     *                 B(1,J),1)
   20    CONTINUE
      ELSE
C
C        Solve A*X = B where A = L*L'.
C
         DO 40 J = 1, NRHS
C
C           Solve L*X = B, overwriting B with X.
C
            CALL ZTBSV('Lower','No transpose','Non-unit',N,KD,AB,LDAB,
     *                 B(1,J),1)
C
C           Solve L'*X = B, overwriting B with X.
C
            CALL ZTBSV('Lower','Conjugate transpose','Non-unit',N,KD,AB,
     *                 LDAB,B(1,J),1)
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F07HSF (ZPBTRS)
C
      END
