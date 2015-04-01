      SUBROUTINE F07QUF(UPLO,N,AP,IPIV,ANORM,RCOND,WORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZSPCON(UPLO,N,AP,IPIV,ANORM,RCOND,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZSPCON estimates the reciprocal of the condition number of a complex
C  symmetric packed matrix A using the factorization A = U*D*U' or
C  A = L*D*L' computed by F07QRF.
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
C  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
C          The block diagonal matrix D and the multipliers used to
C          obtain the factor U or L as computed by F07QRF, stored as a
C          packed triangular matrix.
C
C  IPIV    (input) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D
C          as determined by F07QRF.
C
C  ANORM   (input) REAL
C          The 1-norm of the original matrix A.
C
C  RCOND   (output) REAL
C          The reciprocal of the condition number of the matrix A,
C          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
C          estimate of the 1-norm of inv(A) computed in this routine.
C
C  WORK    (workspace) COMPLEX array, dimension (2*N)
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
      INTEGER           INFO, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), WORK(*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AINVNM
      INTEGER           I, IFAIL, IP, KASE
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F04ZCF, F06AAZ, F07QSF
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
      ELSE IF (ANORM.LT.ZERO) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07QUF/ZSPCON',-INFO)
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
         IP = N*(N+1)/2
         DO 20 I = N, 1, -1
            IF (IPIV(I).GT.0 .AND. AP(IP).EQ.ZERO) RETURN
            IP = IP - I
   20    CONTINUE
      ELSE
C
C        Lower triangular storage: examine D from top to bottom.
C
         IP = 1
         DO 40 I = 1, N
            IF (IPIV(I).GT.0 .AND. AP(IP).EQ.ZERO) RETURN
            IP = IP + N - I + 1
   40    CONTINUE
      END IF
C
C     Estimate the 1-norm of the inverse.
C
      KASE = 0
   60 CONTINUE
      IFAIL = 0
      CALL F04ZCF(KASE,N,WORK,AINVNM,WORK(N+1),IFAIL)
      IF (KASE.NE.0) THEN
C
C        Multiply by inv(L*D*L') or inv(U*D*U').
C
         CALL F07QSF(UPLO,N,1,AP,IPIV,WORK,N,INFO)
         GO TO 60
      END IF
C
C     Compute the estimate of the reciprocal condition number.
C
      IF (AINVNM.NE.ZERO) RCOND = (ONE/AINVNM)/ANORM
C
      RETURN
C
C     End of F07QUF (ZSPCON)
C
      END
