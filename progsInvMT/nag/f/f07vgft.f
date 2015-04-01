      SUBROUTINE F07VGF(NORM,UPLO,DIAG,N,KD,AB,LDAB,RCOND,WORK,IWORK,
     *                  INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DTBCON(NORM,UPLO,DIAG,N,KD,AB,LDAB,RCOND,WORK,
     *                  IWORK,INFO)
C
C  Purpose
C  =======
C
C  DTBCON estimates the reciprocal of the condition number of a
C  triangular band matrix A, in either the 1-norm or the infinity-norm.
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
C  KD      (input) INTEGER
C          The number of superdiagonals or subdiagonals of the
C          triangular band matrix A.  KD >= 0.
C
C  AB      (input) REAL array, dimension (LDAB,N)
C          The upper or lower triangular band matrix A, stored in the
C          first kd+1 rows of the array. The j-th column of A is stored
C          in the j-th column of the array AB as follows:
C          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
C          If DIAG = 'U', the diagonal elements of A are not referenced
C          and are assumed to be 1.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KD+1.
C
C  RCOND   (output) REAL
C          The reciprocal of the condition number of the matrix A,
C          computed as RCOND = 1/(norm(A) * norm(inv(A))).
C
C  WORK    (workspace) REAL array, dimension (3*N)
C
C  IWORK   (workspace) INTEGER array, dimension (N)
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
      INTEGER           INFO, KD, LDAB, N
      CHARACTER         DIAG, NORM, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(LDAB,*), WORK(*)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AINVNM, ANORM, SCALE, SMLNUM, XNORM
      INTEGER           IFAIL, IX, KASE, KASE1
      LOGICAL           NOUNIT, ONENRM, UPPER
      CHARACTER         NORMIN
C     .. External Functions ..
      DOUBLE PRECISION  F06RLF, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          F06RLF, X02AMF, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F04YCF, F06AAZ, F07AGZ, F07VGZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, DBLE
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
      ELSE IF (KD.LT.0) THEN
         INFO = -5
      ELSE IF (LDAB.LT.KD+1) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07VGF/DTBCON',-INFO)
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
      SMLNUM = X02AMF()*DBLE(MAX(1,N))
C
C     Compute the norm of the triangular matrix A.
C
      ANORM = F06RLF(NORM,UPLO,DIAG,N,KD,AB,LDAB,WORK)
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
         CALL F04YCF(KASE,N,WORK,AINVNM,WORK(N+1),IWORK,IFAIL)
         IF (KASE.NE.0) THEN
            IF (KASE.EQ.KASE1) THEN
C
C              Multiply by inv(A).
C
               CALL F07VGZ(UPLO,'No transpose',DIAG,NORMIN,N,KD,AB,LDAB,
     *                     WORK,SCALE,WORK(2*N+1),INFO)
            ELSE
C
C              Multiply by inv(A').
C
               CALL F07VGZ(UPLO,'Transpose',DIAG,NORMIN,N,KD,AB,LDAB,
     *                     WORK,SCALE,WORK(2*N+1),INFO)
            END IF
            NORMIN = 'Y'
C
C           Multiply by 1/SCALE if doing so will not cause overflow.
C
            IF (SCALE.NE.ONE) THEN
               IX = IDAMAX(N,WORK,1)
               XNORM = ABS(WORK(IX))
               IF (SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO) GO TO 40
               CALL F07AGZ(N,SCALE,WORK,1)
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
C     End of F07VGF (DTBCON)
C
      END
