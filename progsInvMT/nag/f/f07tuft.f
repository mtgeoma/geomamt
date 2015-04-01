      SUBROUTINE F07TUF(NORM,UPLO,DIAG,N,A,LDA,RCOND,WORK,RWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZTRCON(NORM,UPLO,DIAG,N,A,LDA,RCOND,WORK,RWORK,
     *                  INFO)
C
C  Purpose
C  =======
C
C  ZTRCON estimates the reciprocal of the condition number of a
C  triangular matrix A, in either the 1-norm or the infinity-norm.
C
C  The norm of A is computed and an estimate is obtained for
C  norm(inv(A)), then the reciprocal of the condition number is
C  computed as
C     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
C
C  Arguments
C  =========
C
C  NORM    (input) CHARACTER*1
C          Specifies whether the 1-norm condition number or the
C          infinity-norm condition number is required:
C          = '1' or 'O':  1-norm
C          = 'I':         Infinity-norm
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the matrix A is upper or lower triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  DIAG    (input) CHARACTER*1
C          Specifies whether or not the matrix A is unit triangular.
C          = 'N':  Non-unit triangular
C          = 'U':  Unit triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
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
C  RCOND   (output) REAL
C          The reciprocal of the condition number of the matrix A,
C          computed as RCOND = 1/(norm(A) * norm(inv(A))).
C
C  WORK    (workspace) COMPLEX array, dimension (2*N)
C
C  RWORK   (workspace) REAL array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -k, the k-th argument had an illegal value
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
      DOUBLE PRECISION  RCOND
      INTEGER           INFO, LDA, N
      CHARACTER         DIAG, NORM, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), WORK(*)
      DOUBLE PRECISION  RWORK(*)
C     .. Local Scalars ..
      COMPLEX*16        ZDUM
      DOUBLE PRECISION  AINVNM, ANORM, SCALE, SMLNUM, XNORM
      INTEGER           IFAIL, IX, KASE, KASE1
      LOGICAL           NOUNIT, ONENRM, UPPER
      CHARACTER         NORMIN
C     .. External Functions ..
      DOUBLE PRECISION  F06UJF, X02ANF
      INTEGER           IZAMAX
      EXTERNAL          F06UJF, X02ANF, IZAMAX
C     .. External Subroutines ..
      EXTERNAL          F04ZCF, F06AAZ, F07AUZ, F07TUZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, MAX, DBLE
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(ZDUM) = ABS(DBLE(ZDUM)) + ABS(DIMAG(ZDUM))
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      ONENRM = NORM .EQ. '1' .OR. (NORM.EQ.'O' .OR. NORM.EQ.'o')
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
      IF ( .NOT. ONENRM .AND. .NOT. (NORM.EQ.'I' .OR. NORM.EQ.'i')) THEN
         INFO = -1
      ELSE IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l'))
     *         THEN
         INFO = -2
      ELSE IF ( .NOT. NOUNIT .AND. .NOT. (DIAG.EQ.'U' .OR. DIAG.EQ.'u'))
     *         THEN
         INFO = -3
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07TUF/ZTRCON',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) THEN
         RCOND = ONE
         RETURN
      END IF
C
      RCOND = ZERO
      SMLNUM = X02ANF()*DBLE(MAX(1,N))
C
C     Compute the norm of the triangular matrix A.
C
      ANORM = F06UJF(NORM,UPLO,DIAG,N,N,A,LDA,RWORK)
C
C     Continue only if ANORM > 0.
C
      IF (ANORM.GT.ZERO) THEN
C
C        Estimate the norm of the inverse of A.
C
         AINVNM = ZERO
         NORMIN = 'N'
         IF (ONENRM) THEN
            KASE1 = 1
         ELSE
            KASE1 = 2
         END IF
         KASE = 0
   20    CONTINUE
         IFAIL = 0
         CALL F04ZCF(KASE,N,WORK,AINVNM,WORK(N+1),IFAIL)
         IF (KASE.NE.0) THEN
            IF (KASE.EQ.KASE1) THEN
C
C              Multiply by inv(A).
C
               CALL F07TUZ(UPLO,'No transpose',DIAG,NORMIN,N,A,LDA,WORK,
     *                     SCALE,RWORK,INFO)
            ELSE
C
C              Multiply by inv(A').
C
               CALL F07TUZ(UPLO,'Conjugate transpose',DIAG,NORMIN,N,A,
     *                     LDA,WORK,SCALE,RWORK,INFO)
            END IF
            NORMIN = 'Y'
C
C           Multiply by 1/SCALE if doing so will not cause overflow.
C
            IF (SCALE.NE.ONE) THEN
               IX = IZAMAX(N,WORK,1)
               XNORM = CABS1(WORK(IX))
               IF (SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO) GO TO 40
               CALL F07AUZ(N,SCALE,WORK,1)
            END IF
            GO TO 20
         END IF
C
C        Compute the estimate of the reciprocal condition number.
C
         IF (AINVNM.NE.ZERO) RCOND = (ONE/ANORM)/AINVNM
      END IF
C
   40 CONTINUE
      RETURN
C
C     End of F07TUF (ZTRCON)
C
      END
