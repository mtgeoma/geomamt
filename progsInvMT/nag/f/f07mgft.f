      SUBROUTINE F07MGF(UPLO,N,A,LDA,IPIV,ANORM,RCOND,WORK,IWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DSYCON(UPLO,N,A,LDA,IPIV,ANORM,RCOND,WORK,IWORK,
     *                  INFO)
C
C  Purpose
C  =======
C
C  DSYCON estimates the reciprocal of the condition number of a real
C  symmetric matrix A using the factorization A = U*D*U' or A = L*D*L'
C  computed by F07MDF.
C
C  An estimate is obtained for norm(inv(A)), and the reciprocal of the
C  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the details of the factorization are stored
C          as an upper or lower triangular matrix.
C          = 'U':  Upper triangular (form is A = U*D*U')
C          = 'L':  Lower triangular (form is A = L*D*L')
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input) REAL array, dimension (LDA,N)
C          The block diagonal matrix D and the multipliers used to
C          obtain the factor U or L as computed by F07MDF.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (input) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D
C          as determined by F07MDF.
C
C  ANORM   (input) REAL
C          The 1-norm of the original matrix A.
C
C  RCOND   (output) REAL
C          The reciprocal of the condition number of the matrix A,
C          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
C          estimate of the 1-norm of inv(A) computed in this routine.
C
C  WORK    (workspace) REAL array, dimension (2*N)
C
C  IWORK    (workspace) INTEGER array, dimension (N)
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
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANORM, RCOND
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(*)
      INTEGER           IPIV(*), IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AINVNM
      INTEGER           I, IFAIL, KASE
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F04YCF, F06AAZ, F07MEF
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
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (ANORM.LT.ZERO) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07MGF/DSYCON',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      RCOND = ZERO
      IF (N.EQ.0) THEN
         RCOND = ONE
         RETURN
      ELSE IF (ANORM.LE.ZERO) THEN
         RETURN
      END IF
C
C     Check that the diagonal matrix D is nonsingular.
C
      IF (UPPER) THEN
C
C        Upper triangular storage: examine D from bottom to top
C
         DO 20 I = N, 1, -1
            IF (IPIV(I).GT.0 .AND. A(I,I).EQ.ZERO) RETURN
   20    CONTINUE
      ELSE
C
C        Lower triangular storage: examine D from top to bottom.
C
         DO 40 I = 1, N
            IF (IPIV(I).GT.0 .AND. A(I,I).EQ.ZERO) RETURN
   40    CONTINUE
      END IF
C
C     Estimate the 1-norm of the inverse.
C
      KASE = 0
   60 CONTINUE
      IFAIL = 0
      CALL F04YCF(KASE,N,WORK,AINVNM,WORK(N+1),IWORK,IFAIL)
      IF (KASE.NE.0) THEN
C
C        Multiply by inv(L*D*L') or inv(U*D*U').
C
         CALL F07MEF(UPLO,N,1,A,LDA,IPIV,WORK,N,INFO)
         GO TO 60
      END IF
C
C     Compute the estimate of the reciprocal condition number.
C
      IF (AINVNM.NE.ZERO) RCOND = (ONE/AINVNM)/ANORM
C
      RETURN
C
C     End of F07MGF (DSYCON)
C
      END
