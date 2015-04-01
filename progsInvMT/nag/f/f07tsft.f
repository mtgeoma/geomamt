      SUBROUTINE F07TSF(UPLO,TRANS,DIAG,N,NRHS,A,LDA,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1025 (JUN 1993).
C     .. Entry Points ..
      ENTRY             ZTRTRS(UPLO,TRANS,DIAG,N,NRHS,A,LDA,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  ZTRTRS solves a complex triangular system of linear equations with
C  multiple right hand sides, of the form
C
C     A * X = B,  A**T * X = B,  or  A**H * X = B,
C
C  where A is a triangular matrix of order N, A**T is the transpose of
C  A, A**H is the conjugate transpose of A, and B is an N by NRHS
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
C          Specifies the form of the equations.
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
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrix B.  NRHS >= 0.
C
C  A       (input) COMPLEX array, dimension (LDA,N)
C          The triangular matrix A.  If UPLO = 'U', the leading n by n
C          upper triangular part of the array A contains the upper
C          triangular matrix, and the strictly lower triangular part of
C          A is not referenced.  If UPLO = 'L', the leading n by n lower
C          triangular part of the array A contains the lower triangular
C          matrix, and the strictly upper triangular part of A is not
C          referenced.  If DIAG = 'U', the diagonal elements of A are
C          also not referenced and are assumed to be 1.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
C          On entry, the right hand side matrix B.
C          On exit, if INFO = 0, the solution matrix X.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the k-th diagonal element of A is zero,
C               indicating that the matrix is singular and the solution
C               X has not been computed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, N, NRHS
      CHARACTER         DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*)
C     .. Local Scalars ..
      LOGICAL           NOUNIT
C     .. External Subroutines ..
      EXTERNAL          ZTRSM, ZTRSV, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
      IF ( .NOT. (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
     *    .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
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
      ELSE IF (NRHS.LT.0) THEN
         INFO = -5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -7
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -9
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07TSF/ZTRTRS',-INFO)
         RETURN
      END IF
C
C     Check for singularity.
C
      IF (NOUNIT) THEN
         DO 20 INFO = 1, N
            IF (A(INFO,INFO).EQ.ZERO) RETURN
   20    CONTINUE
      END IF
      INFO = 0
C
C     Solve A * X = B,  A**T * X = B,  or  A**H * X = B.
C
      IF (NRHS.EQ.1) THEN
         CALL ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,B,1)
      ELSE
         CALL ZTRSM('Left',UPLO,TRANS,DIAG,N,NRHS,ONE,A,LDA,B,LDB)
      END IF
C
      RETURN
C
C     End of F07TSF (ZTRTRS)
C
      END
