      SUBROUTINE F07VSF(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZTBTRS(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,B,LDB,
     *                  INFO)
C
C  Purpose
C  =======
C
C  ZTBTRS solves a triangular system of the form
C
C     A * X = B,  A**T * X = B,  or  A**H * X = B,
C
C  where A is a triangular band matrix of order N, A**T is the transpose
C  of A, A**H is the conjugate transpose of A, and B is an N by NRHS
C  matrix.  A check is made to verify that A is nonsingular.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the matrix A is upper or lower triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  TRANS   (input) CHARACTER*1
C          Specifies the operation applied to A.
C          = 'N':  Solve  A * X = B     (No transpose)
C          = 'T':  Solve  A**T * X = B  (Transpose)
C          = 'C':  Solve  A**H * X = B  (Conjugate transpose)
C
C  DIAG    (input) CHARACTER*1
C          Specifies whether or not the matrix A is unit triangular.
C          = 'N':  Non-unit triangular
C          = 'U':  Unit triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  KD      (input) INTEGER
C          The number of superdiagonals or subdiagonals of the
C          triangular band matrix A.  KD >= 0.
C
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrix B.  NRHS >= 0.
C
C  AB      (input) COMPLEX array, dimension (LDAB,N)
C          The upper or lower triangular band matrix A, stored in the
C          first kd+1 rows of AB.  The j-th column of A is stored
C          in the j-th column of the array AB as follows:
C          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
C          If DIAG = 'U', the diagonal elements of A are not referenced
C          and are assumed to be 1.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KD+1.
C
C  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
C          On entry, the right hand side vectors B for the system of
C          linear equations.
C          On exit, if INFO = 0, the solution vectors X.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the k-th diagonal element of A is zero,
C               indicating that the matrix is singular and the solutions
C               X have not been computed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, KD, LDAB, LDB, N, NRHS
      CHARACTER         DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX*16        AB(LDAB,*), B(LDB,*)
C     .. Local Scalars ..
      INTEGER           J
      LOGICAL           NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          ZTBSV, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF ( .NOT. (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
     *         .AND. .NOT. (TRANS.EQ.'T' .OR. TRANS.EQ.'t')
     *         .AND. .NOT. (TRANS.EQ.'C' .OR. TRANS.EQ.'c')) THEN
         INFO = -2
      ELSE IF ( .NOT. NOUNIT .AND. .NOT. (DIAG.EQ.'U' .OR. DIAG.EQ.'u'))
     *         THEN
         INFO = -3
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (KD.LT.0) THEN
         INFO = -5
      ELSE IF (NRHS.LT.0) THEN
         INFO = -6
      ELSE IF (LDAB.LT.KD+1) THEN
         INFO = -8
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -10
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07VSF/ZTBTRS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Check for singularity.
C
      IF (NOUNIT) THEN
         IF (UPPER) THEN
            DO 20 INFO = 1, N
               IF (AB(KD+1,INFO).EQ.ZERO) RETURN
   20       CONTINUE
         ELSE
            DO 40 INFO = 1, N
               IF (AB(1,INFO).EQ.ZERO) RETURN
   40       CONTINUE
         END IF
      END IF
      INFO = 0
C
C     Solve A * X = B,  A**T * X = B,  or  A**H * X = B.
C
      DO 60 J = 1, NRHS
         CALL ZTBSV(UPLO,TRANS,DIAG,N,KD,AB,LDAB,B(1,J),1)
   60 CONTINUE
C
      RETURN
C
C     End of F07VSF (ZTBTRS)
C
      END
