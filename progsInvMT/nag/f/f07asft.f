      SUBROUTINE F07ASF(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1019 (JUN 1993).
C     .. Entry Points ..
      ENTRY             ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  ZGETRS solves a system of linear equations
C     A * X = B,  A**T * X = B,  or  A**H * X = B
C  with a general n by n matrix A using the LU factorization computed
C  by F07ARF.
C
C  Arguments
C  =========
C
C  TRANS   (input) CHARACTER*1
C          Specifies the form of the system of equations.
C          = 'N':  A * X = B     (No transpose)
C          = 'T':  A**T * X = B  (Transpose)
C          = 'C':  A**H * X = B  (Conjugate transpose)
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrix B.  NRHS >= 0.
C
C  A       (input) COMPLEX array, dimension (LDA,N)
C          The factors L and U from the factorization A = P*L*U
C          as computed by F07ARF.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (input) INTEGER array, dimension (N)
C          The pivot indices from F07ARF; for 1<=i<=N, row i of the
C          matrix was interchanged with row IPIV(i).
C
C  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
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
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, N, NRHS
      CHARACTER         TRANS
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      LOGICAL           NOTRAN
C     .. External Subroutines ..
      EXTERNAL          ZTRSM, ZTRSV, F06AAZ, F07ARY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
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
      ELSE IF (NRHS.LT.0) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07ASF/ZGETRS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
C
      IF (NOTRAN) THEN
C
C        Solve A * X = B.
C
C        Apply row interchanges to the right hand sides.
C
         CALL F07ARY(NRHS,B,LDB,1,N,IPIV,1)
C
C        Solve L*Y = B, overwriting B with Y, and then solve U*X = Y,
C        overwriting Y with X.
C
         IF (NRHS.EQ.1) THEN
            CALL ZTRSV('Lower','No transpose','Unit',N,A,LDA,B,1)
            CALL ZTRSV('Upper','No transpose','Non-unit',N,A,LDA,B,1)
         ELSE
            CALL ZTRSM('Left','Lower','No transpose','Unit',N,NRHS,ONE,
     *                 A,LDA,B,LDB)
            CALL ZTRSM('Left','Upper','No transpose','Non-unit',N,NRHS,
     *                 ONE,A,LDA,B,LDB)
         END IF
      ELSE
C
C        Solve A**T * X = B  or A**H * X = B.
C
C        Solve U'*Y = B, overwriting B with Y, and then solve L'*X = Y,
C        overwriting Y with X.
C
         IF (NRHS.EQ.1) THEN
            CALL ZTRSV('Upper',TRANS,'Non-unit',N,A,LDA,B,1)
            CALL ZTRSV('Lower',TRANS,'Unit',N,A,LDA,B,1)
         ELSE
            CALL ZTRSM('Left','Upper',TRANS,'Non-unit',N,NRHS,ONE,A,LDA,
     *                 B,LDB)
            CALL ZTRSM('Left','Lower',TRANS,'Unit',N,NRHS,ONE,A,LDA,B,
     *                 LDB)
         END IF
C
C        Apply row interchanges to the solution vectors.
C
         CALL F07ARY(NRHS,B,LDB,1,N,IPIV,-1)
      END IF
C
      RETURN
C
C     End of F07ASF (ZGETRS)
C
      END
