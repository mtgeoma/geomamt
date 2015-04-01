      SUBROUTINE F07BEF(TRANS,N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DGBTRS(TRANS,N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,
     *                  INFO)
C
C  Purpose
C  =======
C
C  DGBTRS solves a system of linear equations
C     A * X = B  or  A' * X = B
C  with a general band matrix A using the LU factorization computed
C  by F07BDF.
C
C  Arguments
C  =========
C
C  TRANS   (input) CHARACTER*1
C          Specifies the form of the system of equations.
C          = 'N':  A * X = B  (No transpose)
C          = 'T':  A'* X = B  (Transpose)
C          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
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
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrix B.  NRHS >= 0.
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
C  B       (input/output) REAL array, dimension (LDB,NRHS)
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
      INTEGER           INFO, KL, KU, LDAB, LDB, N, NRHS
      CHARACTER         TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(LDAB,*), B(LDB,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           I, J, KD, L, LM
      LOGICAL           LNOTI, NOTRAN
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, DAXPY, DTBSV
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      NOTRAN = (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
      IF ( .NOT. NOTRAN .AND. .NOT. (TRANS.EQ.'T' .OR. TRANS.EQ.'t')
     *     .AND. .NOT. (TRANS.EQ.'C' .OR. TRANS.EQ.'c')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (KL.LT.0) THEN
         INFO = -3
      ELSE IF (KU.LT.0) THEN
         INFO = -4
      ELSE IF (NRHS.LT.0) THEN
         INFO = -5
      ELSE IF (LDAB.LT.(2*KL+KU+1)) THEN
         INFO = -7
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -10
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07BEF/DGBTRS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
C
      KD = KU + KL + 1
      LNOTI = KL .GT. 0
C
      IF (NOTRAN) THEN
C
C        Solve  A*X = B.
C
         DO 40 I = 1, NRHS
C
C           Solve L*X = B, overwriting B with X.
C
C           L is represented as a product of permutations and unit lower
C           triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
C           where each transformation L(i) is a rank-one modification of
C           the identity matrix.
C
            IF (LNOTI) THEN
               DO 20 J = 1, N - 1
                  LM = MIN(KL,N-J)
                  L = IPIV(J)
                  T = B(L,I)
                  IF (L.NE.J) THEN
                     B(L,I) = B(J,I)
                     B(J,I) = T
                  END IF
                  CALL DAXPY(LM,-T,AB(KD+1,J),1,B(J+1,I),1)
   20          CONTINUE
            END IF
C
C           Solve U*X = B, overwriting B with X.
C
            CALL DTBSV('Upper','No transpose','Non-unit',N,KL+KU,AB,
     *                 LDAB,B(1,I),1)
   40    CONTINUE
C
      ELSE
C
C        Solve A'*X = B.
C
         DO 80 I = 1, NRHS
C
C           Solve U'*X = B, overwriting B with X.
C
            CALL DTBSV('Upper','Transpose','Non-unit',N,KL+KU,AB,LDAB,
     *                 B(1,I),1)
C
C           Solve L'*X = B, overwriting B with X.
C
            IF (LNOTI) THEN
               DO 60 J = N - 1, 1, -1
                  LM = MIN(KL,N-J)
                  B(J,I) = B(J,I) - DDOT(LM,AB(KD+1,J),1,B(J+1,I),1)
                  L = IPIV(J)
                  IF (L.NE.J) THEN
                     T = B(L,I)
                     B(L,I) = B(J,I)
                     B(J,I) = T
                  END IF
   60          CONTINUE
            END IF
   80    CONTINUE
      END IF
      RETURN
C
C     End of F07BEF (DGBTRS)
C
      END
