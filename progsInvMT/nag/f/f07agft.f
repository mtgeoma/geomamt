      SUBROUTINE F07AGF(NORM,N,A,LDA,ANORM,RCOND,WORK,IWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DGECON(NORM,N,A,LDA,ANORM,RCOND,WORK,IWORK,INFO)
C
C  Purpose
C  =======
C
C  DGECON estimates the reciprocal of the condition number of a real
C  general matrix A, in either the 1-norm or the infinity-norm, using
C  the LU factorization computed by F07ADF.
C
C  An estimate is obtained for norm(inv(A)), and the reciprocal of the
C  condition number is computed as
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
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input) REAL array, dimension (LDA,N)
C          The factors L and U from the factorization A = P*L*U
C          as computed by F07ADF.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  ANORM   (input) REAL
C          If NORM = '1' or 'O', the 1-norm of the original matrix A.
C          If NORM = 'I', the infinity-norm of the original matrix A.
C
C  RCOND   (output) REAL
C          The reciprocal of the condition number of the matrix A,
C          computed as RCOND = 1/(norm(A) * norm(inv(A))).
C
C  WORK    (workspace) REAL array, dimension (4*N)
C
C  IWORK   (workspace) INTEGER array, dimension (N)
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
      CHARACTER         NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(*)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AINVNM, SCALE, SL, SMLNUM, SU
      INTEGER           IFAIL, IX, KASE, KASE1
      LOGICAL           ONENRM
      CHARACTER         NORMIN
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           IDAMAX
      EXTERNAL          X02AMF, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F04YCF, F06AAZ, F07AGZ, F07TGZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      ONENRM = NORM .EQ. '1' .OR. (NORM.EQ.'O' .OR. NORM.EQ.'o')
      IF ( .NOT. ONENRM .AND. .NOT. (NORM.EQ.'I' .OR. NORM.EQ.'i')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (ANORM.LT.ZERO) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07AGF/DGECON',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      RCOND = ZERO
      IF (N.EQ.0) THEN
         RCOND = ONE
         RETURN
      ELSE IF (ANORM.EQ.ZERO) THEN
         RETURN
      END IF
C
      SMLNUM = X02AMF()
C
C     Estimate the norm of inv(A).
C
      AINVNM = ZERO
      NORMIN = 'N'
      IF (ONENRM) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   20 CONTINUE
      IFAIL = 0
      CALL F04YCF(KASE,N,WORK,AINVNM,WORK(N+1),IWORK,IFAIL)
      IF (KASE.NE.0) THEN
         IF (KASE.EQ.KASE1) THEN
C
C           Multiply by inv(L).
C
            CALL F07TGZ('Lower','No transpose','Unit',NORMIN,N,A,LDA,
     *                  WORK,SL,WORK(2*N+1),INFO)
C
C           Multiply by inv(U).
C
            CALL F07TGZ('Upper','No transpose','Non-unit',NORMIN,N,A,
     *                  LDA,WORK,SU,WORK(3*N+1),INFO)
         ELSE
C
C           Multiply by inv(U').
C
            CALL F07TGZ('Upper','Transpose','Non-unit',NORMIN,N,A,LDA,
     *                  WORK,SU,WORK(3*N+1),INFO)
C
C           Multiply by inv(L').
C
            CALL F07TGZ('Lower','Transpose','Unit',NORMIN,N,A,LDA,WORK,
     *                  SL,WORK(2*N+1),INFO)
         END IF
C
C        Divide X by 1/(SL*SU) if doing so will not cause overflow.
C
         SCALE = SL*SU
         NORMIN = 'Y'
         IF (SCALE.NE.ONE) THEN
            IX = IDAMAX(N,WORK,1)
            IF (SCALE.LT.ABS(WORK(IX))*SMLNUM .OR. SCALE.EQ.ZERO)
     *          GO TO 40
            CALL F07AGZ(N,SCALE,WORK,1)
         END IF
         GO TO 20
      END IF
C
C     Compute the estimate of the reciprocal condition number.
C
      IF (AINVNM.NE.ZERO) RCOND = (ONE/AINVNM)/ANORM
C
   40 CONTINUE
      RETURN
C
C     End of F07AGF (DGECON)
C
      END
