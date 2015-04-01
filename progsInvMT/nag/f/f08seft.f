      SUBROUTINE F08SEF(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DSYGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  DSYGST reduces a real symmetric-definite generalized eigenproblem
C  to standard form.
C
C  If ITYPE = 1, the problem is A*x = lambda*B*x,
C  and A is overwritten by inv(U')*A*inv(U) or inv(L)*A*inv(L')
C
C  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
C  B*A*x = lambda*x, and A is overwritten by U*A*U` or L'*A*L.
C
C  B must have been previously factorized as U'*U or L*L' by DPOTRF.
C
C  Arguments
C  =========
C
C  ITYPE   (input) INTEGER
C          = 1: compute inv(U')*A*inv(U) or inv(L)*A*inv(L');
C          = 2 or 3: compute U*A*U' or L'*A*L.
C
C  UPLO    (input) CHARACTER
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored, and how B has been factorized.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrices A and B.  N >= 0.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n by n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n by n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C
C          On exit, if INFO = 0, the transformed matrix, stored in the
C          same format as A.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
C          The triangular factor from the Cholesky factorization of B,
C          as returned by DPOTRF.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit.
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, HALF
      PARAMETER         (ONE=1.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, ITYPE, LDA, LDB, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*)
C     .. Local Scalars ..
      INTEGER           K, KB, NB
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          DSYMM, DSYR2K, DTRMM, DTRSM, F06AAZ, F07ZAZ,
     *                  F08SEZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF (ITYPE.LT.1 .OR. ITYPE.GT.3) THEN
         INFO = -1
      ELSE IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l'))
     *         THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08SEF/DSYGST',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F08SEF',NB,0)
C
      IF (NB.LE.1 .OR. NB.GE.N) THEN
C
C        Use unblocked code
C
         CALL F08SEZ(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
      ELSE
C
C        Use blocked code
C
         IF (ITYPE.EQ.1) THEN
            IF (UPPER) THEN
C
C              Compute inv(U')*A*inv(U)
C
               DO 20 K = 1, N, NB
                  KB = MIN(N-K+1,NB)
C
C                 Update the upper triangle of A(k:n,k:n)
C
                  CALL F08SEZ(ITYPE,UPLO,KB,A(K,K),LDA,B(K,K),LDB,INFO)
                  IF (K+KB.LE.N) THEN
                     CALL DTRSM('Left',UPLO,'Transpose','Non-unit',KB,
     *                          N-K-KB+1,ONE,B(K,K),LDB,A(K,K+KB),LDA)
                     CALL DSYMM('Left',UPLO,KB,N-K-KB+1,-HALF,A(K,K),
     *                          LDA,B(K,K+KB),LDB,ONE,A(K,K+KB),LDA)
                     CALL DSYR2K(UPLO,'Transpose',N-K-KB+1,KB,-ONE,
     *                           A(K,K+KB),LDA,B(K,K+KB),LDB,ONE,
     *                           A(K+KB,K+KB),LDA)
                     CALL DSYMM('Left',UPLO,KB,N-K-KB+1,-HALF,A(K,K),
     *                          LDA,B(K,K+KB),LDB,ONE,A(K,K+KB),LDA)
                     CALL DTRSM('Right',UPLO,'No transpose','Non-unit',
     *                          KB,N-K-KB+1,ONE,B(K+KB,K+KB),LDB,
     *                          A(K,K+KB),LDA)
                  END IF
   20          CONTINUE
            ELSE
C
C              Compute inv(L)*A*inv(L')
C
               DO 40 K = 1, N, NB
                  KB = MIN(N-K+1,NB)
C
C                 Update the lower triangle of A(k:n,k:n)
C
                  CALL F08SEZ(ITYPE,UPLO,KB,A(K,K),LDA,B(K,K),LDB,INFO)
                  IF (K+KB.LE.N) THEN
                     CALL DTRSM('Right',UPLO,'Transpose','Non-unit',
     *                          N-K-KB+1,KB,ONE,B(K,K),LDB,A(K+KB,K),
     *                          LDA)
                     CALL DSYMM('Right',UPLO,N-K-KB+1,KB,-HALF,A(K,K),
     *                          LDA,B(K+KB,K),LDB,ONE,A(K+KB,K),LDA)
                     CALL DSYR2K(UPLO,'No transpose',N-K-KB+1,KB,-ONE,
     *                           A(K+KB,K),LDA,B(K+KB,K),LDB,ONE,
     *                           A(K+KB,K+KB),LDA)
                     CALL DSYMM('Right',UPLO,N-K-KB+1,KB,-HALF,A(K,K),
     *                          LDA,B(K+KB,K),LDB,ONE,A(K+KB,K),LDA)
                     CALL DTRSM('Left',UPLO,'No transpose','Non-unit',
     *                          N-K-KB+1,KB,ONE,B(K+KB,K+KB),LDB,
     *                          A(K+KB,K),LDA)
                  END IF
   40          CONTINUE
            END IF
         ELSE
            IF (UPPER) THEN
C
C              Compute U*A*U'
C
               DO 60 K = 1, N, NB
                  KB = MIN(N-K+1,NB)
C
C                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
C
                  CALL DTRMM('Left',UPLO,'No transpose','Non-unit',K-1,
     *                       KB,ONE,B,LDB,A(1,K),LDA)
                  CALL DSYMM('Right',UPLO,K-1,KB,HALF,A(K,K),LDA,B(1,K),
     *                       LDB,ONE,A(1,K),LDA)
                  CALL DSYR2K(UPLO,'No transpose',K-1,KB,ONE,A(1,K),LDA,
     *                        B(1,K),LDB,ONE,A,LDA)
                  CALL DSYMM('Right',UPLO,K-1,KB,HALF,A(K,K),LDA,B(1,K),
     *                       LDB,ONE,A(1,K),LDA)
                  CALL DTRMM('Right',UPLO,'Transpose','Non-unit',K-1,KB,
     *                       ONE,B(K,K),LDB,A(1,K),LDA)
                  CALL F08SEZ(ITYPE,UPLO,KB,A(K,K),LDA,B(K,K),LDB,INFO)
   60          CONTINUE
            ELSE
C
C              Compute L'*A*L
C
               DO 80 K = 1, N, NB
                  KB = MIN(N-K+1,NB)
C
C                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
C
                  CALL DTRMM('Right',UPLO,'No transpose','Non-unit',KB,
     *                       K-1,ONE,B,LDB,A(K,1),LDA)
                  CALL DSYMM('Left',UPLO,KB,K-1,HALF,A(K,K),LDA,B(K,1),
     *                       LDB,ONE,A(K,1),LDA)
                  CALL DSYR2K(UPLO,'Transpose',K-1,KB,ONE,A(K,1),LDA,
     *                        B(K,1),LDB,ONE,A,LDA)
                  CALL DSYMM('Left',UPLO,KB,K-1,HALF,A(K,K),LDA,B(K,1),
     *                       LDB,ONE,A(K,1),LDA)
                  CALL DTRMM('Left',UPLO,'Transpose','Non-unit',KB,K-1,
     *                       ONE,B(K,K),LDB,A(K,1),LDA)
                  CALL F08SEZ(ITYPE,UPLO,KB,A(K,K),LDA,B(K,K),LDB,INFO)
   80          CONTINUE
            END IF
         END IF
      END IF
      RETURN
C
C     End of F08SEF (DSYGST)
C
      END
