      SUBROUTINE F07USF(UPLO,TRANS,DIAG,N,NRHS,AP,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZTPTRS(UPLO,TRANS,DIAG,N,NRHS,AP,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  ZTPTRS solves a complex triangular system of linear equations with
C  multiple right hand sides, of the form
C
C     A * X = B,  A**T * X = B,  or  A**H * X = B,
C
C  where A is a triangular matrix of order N stored in packed format,
C  A**T is the transpose of A, A**H is the conjugate transpose of A,
C  and B is an N by NRHS matrix.  A check is made to verify that A is
C  nonsingular.
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
C  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
C          The upper or lower triangular matrix A, packed columnwise in
C          a linear array:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i.
C          If DIAG = 'U', the diagonal elements of A are not referenced
C          and are assumed to be 1.
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
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDB, N, NRHS
      CHARACTER         DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), B(LDB,*)
C     .. Local Scalars ..
      INTEGER           J, JC
      LOGICAL           NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          ZTPSV, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
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
      ELSE IF (NRHS.LT.0) THEN
         INFO = -5
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07USF/ZTPTRS',-INFO)
         RETURN
      END IF
C
C     Check for singularity.
C
      IF (NOUNIT) THEN
         IF (UPPER) THEN
            JC = 1
            DO 20 INFO = 1, N
               IF (AP(JC+INFO-1).EQ.ZERO) RETURN
               JC = JC + INFO
   20       CONTINUE
         ELSE
            JC = 1
            DO 40 INFO = 1, N
               IF (AP(JC).EQ.ZERO) RETURN
               JC = JC + N - INFO + 1
   40       CONTINUE
         END IF
      END IF
      INFO = 0
C
C     Solve  A * X = B,  A**T * X = B,  or  A**H * X = B.
C
      DO 60 J = 1, NRHS
         CALL ZTPSV(UPLO,TRANS,DIAG,N,AP,B(1,J),1)
   60 CONTINUE
C
      RETURN
C
C     End of F07USF (ZTPTRS)
C
      END
