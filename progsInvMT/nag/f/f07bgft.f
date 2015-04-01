      SUBROUTINE F07BGF(NORM,N,KL,KU,AB,LDAB,IPIV,ANORM,RCOND,WORK,
     *                  IWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DGBCON(NORM,N,KL,KU,AB,LDAB,IPIV,ANORM,RCOND,
     *                  WORK,IWORK,INFO)
C
C  Purpose
C  =======
C
C  DGBCON estimates the reciprocal of the condition number of a real
C  general band matrix A, in either the 1-norm or the infinity-norm,
C  using the LU factorization computed by F07BDF.
C
C  An estimate is obtained for norm(inv(A)), and RCOND is computed as
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
C  KL      (input) INTEGER
C          The number of subdiagonals within the band of A.  KL >= 0.
C
C  KU      (input) INTEGER
C          The number of superdiagonals within the band of A.  KU >= 0.
C
C  AB      (input) REAL array, dimension (LDAB,N)
C          Details of the LU factorization of the band matrix A, as
C          computed by F07BDF.  U is stored as an upper triangular band
C          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
C          the multipliers used during the factorization are stored in
C          rows KL+KU+2 to 2*KL+KU+1.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
C
C  IPIV    (input) INTEGER array, dimension (N)
C          The pivot indices; for 1 <= i <= N, row i of the matrix was
C          interchanged with row IPIV(i).
C
C  ANORM   (input) REAL
C          If NORM = '1' or 'O', the 1-norm of the original matrix A.
C          If NORM = 'I', the infinity-norm of the original matrix A.
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
      INTEGER           INFO, KL, KU, LDAB, N
      CHARACTER         NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(LDAB,*), WORK(*)
      INTEGER           IPIV(*), IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AINVNM, SCALE, SMLNUM, T
      INTEGER           IFAIL, IX, J, JP, KASE, KASE1, KD, LM
      LOGICAL           LNOTI, ONENRM
      CHARACTER         NORMIN
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          DDOT, X02AMF, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F04YCF, F06AAZ, F07AGZ, F07VGZ, DAXPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
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
      ELSE IF (KL.LT.0) THEN
         INFO = -3
      ELSE IF (KU.LT.0) THEN
         INFO = -4
      ELSE IF (LDAB.LT.2*KL+KU+1) THEN
         INFO = -6
      ELSE IF (ANORM.LT.ZERO) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07BGF/DGBCON',-INFO)
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
      KD = KL + KU + 1
      LNOTI = KL .GT. 0
      KASE = 0
   20 CONTINUE
      IFAIL = 0
      CALL F04YCF(KASE,N,WORK,AINVNM,WORK(N+1),IWORK,IFAIL)
      IF (KASE.NE.0) THEN
         IF (KASE.EQ.KASE1) THEN
C
C           Multiply by inv(L).
C
            IF (LNOTI) THEN
               DO 40 J = 1, N - 1
                  LM = MIN(KL,N-J)
                  JP = IPIV(J)
                  T = WORK(JP)
                  IF (JP.NE.J) THEN
                     WORK(JP) = WORK(J)
                     WORK(J) = T
                  END IF
                  CALL DAXPY(LM,-T,AB(KD+1,J),1,WORK(J+1),1)
   40          CONTINUE
            END IF
C
C           Multiply by inv(U).
C
            CALL F07VGZ('Upper','No transpose','Non-unit',NORMIN,N,
     *                  KL+KU,AB,LDAB,WORK,SCALE,WORK(2*N+1),INFO)
         ELSE
C
C           Multiply by inv(U').
C
            CALL F07VGZ('Upper','Transpose','Non-unit',NORMIN,N,KL+KU,
     *                  AB,LDAB,WORK,SCALE,WORK(2*N+1),INFO)
C
C           Multiply by inv(L').
C
            IF (LNOTI) THEN
               DO 60 J = N - 1, 1, -1
                  LM = MIN(KL,N-J)
                  WORK(J) = WORK(J) - DDOT(LM,AB(KD+1,J),1,WORK(J+1),1)
                  JP = IPIV(J)
                  IF (JP.NE.J) THEN
                     T = WORK(JP)
                     WORK(JP) = WORK(J)
                     WORK(J) = T
                  END IF
   60          CONTINUE
            END IF
         END IF
C
C        Divide X by 1/SCALE if doing so will not cause overflow.
C
         NORMIN = 'Y'
         IF (SCALE.NE.ONE) THEN
            IX = IDAMAX(N,WORK,1)
            IF (SCALE.LT.ABS(WORK(IX))*SMLNUM .OR. SCALE.EQ.ZERO)
     *          GO TO 80
            CALL F07AGZ(N,SCALE,WORK,1)
         END IF
         GO TO 20
      END IF
C
C     Compute the estimate of the reciprocal condition number.
C
      IF (AINVNM.NE.ZERO) RCOND = (ONE/AINVNM)/ANORM
C
   80 CONTINUE
      RETURN
C
C     End of F07BGF (DGBCON)
C
      END
