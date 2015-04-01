      SUBROUTINE F07MHF(UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX,FERR,
     *                  BERR,WORK,IWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DSYRFS(UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,
     *                  LDX,FERR,BERR,WORK,IWORK,INFO)
C
C  Purpose
C  =======
C
C  DSYRFS improves the computed solution to a system of linear
C  equations when the coefficient matrix is symmetric indefinite and
C  provides error bounds and backward error estimates for the solutions.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrices B and X.  NRHS >= 0.
C
C  A       (input) REAL array, dimension (LDA,N)
C          The symmetric matrix A.  If UPLO = 'U', the leading n by n
C          upper triangular part of A contains the upper triangular part
C          of the matrix A, and the strictly lower triangular part of A
C          is not referenced.  If UPLO = 'L', the leading n by n lower
C          triangular part of A contains the lower triangular part of
C          the matrix A, and the strictly upper triangular part of A is
C          not referenced.
C
C  LDA     (input) INTEGER.
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  AF      (input) REAL array, dimension (LDAF,N)
C          The factored form of the matrix A.  AF contains the block
C          diagonal matrix D and the multipliers used to obtain the
C          factor U or L from the factorization A = U*D*U' or A = L*D*L'
C          as computed by F07MDF.
C
C  LDAF    (input) INTEGER
C          The leading dimension of the array AF.  LDAF >= max(1,N).
C
C  IPIV    (input) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D
C          as determined by F07MDF.
C
C  B       (input) REAL array, dimension (LDB,NRHS)
C          The right hand side vectors for the system of linear
C          equations.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  X       (input/output) REAL array, dimension (LDX,NRHS)
C          On entry, the solution vectors.
C          On exit, the improved solution vectors.
C
C  LDX     (input) INTEGER
C          The leading dimension of the array X.  LDX >= max(1,N).
C
C  FERR    (output) REAL array, dimension (NRHS)
C          The estimated forward error bounds for each solution vector
C          X.  If XTRUE is the true solution, FERR bounds the magnitude
C          of the largest entry in (X - XTRUE) divided by the magnitude
C          of the largest entry in X.  The quality of the error bound
C          depends on the quality of the estimate of norm(inv(A))
C          computed in the code; if the estimate of norm(inv(A)) is
C          accurate, the error bound is guaranteed.
C
C  BERR    (output) REAL array, dimension (NRHS)
C          The componentwise relative backward error of each solution
C          vector (i.e., the smallest relative change in any entry of A
C          or B that makes X an exact solution).
C
C  WORK    (workspace) REAL array, dimension (3*N)
C
C  IWORK   (workspace) INTEGER array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C
C  Further Details
C  ======= =======
C
C  ITMAX is the maximum number of steps of iterative refinement.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      INTEGER           ITMAX
      PARAMETER         (ITMAX=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D+0)
      DOUBLE PRECISION  THREE
      PARAMETER         (THREE=3.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDAF, LDB, LDX, N, NRHS
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AF(LDAF,*), B(LDB,*), BERR(*),
     *                  FERR(*), WORK(*), X(LDX,*)
      INTEGER           IPIV(*), IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, LSTRES, S, XK
      INTEGER           COUNT, I, IFAIL, J, K, KASE, NZ
      LOGICAL           UPPER
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F04YCF, F06AAZ, F07MEF, DAXPY, DCOPY, DSYMV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
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
      ELSE IF (NRHS.LT.0) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDAF.LT.MAX(1,N)) THEN
         INFO = -7
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -10
      ELSE IF (LDX.LT.MAX(1,N)) THEN
         INFO = -12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07MHF/DSYRFS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. NRHS.EQ.0) THEN
         DO 20 J = 1, NRHS
            FERR(J) = ZERO
            BERR(J) = ZERO
   20    CONTINUE
         RETURN
      END IF
C
      EPS = X02AJF()
C
C     NZ = maximum number of nonzero entries in each row of A, plus 1
C
      NZ = N + 1
C
C     Do for each right hand side
C
      DO 280 J = 1, NRHS
C
         COUNT = 1
         LSTRES = THREE
   40    CONTINUE
C
C        Loop until stopping criterion is satisfied.
C
C        Compute residual R = B - A * X
C
         CALL DCOPY(N,B(1,J),1,WORK(N+1),1)
         CALL DSYMV(UPLO,N,-ONE,A,LDA,X(1,J),1,ONE,WORK(N+1),1)
C
C        Compute componentwise relative backward error from formula
C
C        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )
C
C        where 0/0 is treated as 0, and abs(Z) is the componentwise
C        absolute value of the matrix or vector Z.
C
         DO 60 I = 1, N
            WORK(I) = ABS(B(I,J))
   60    CONTINUE
C
C        Compute abs(A)*abs(X) + abs(B).
C
         IF (UPPER) THEN
            DO 100 K = 1, N
               S = ZERO
               XK = ABS(X(K,J))
               DO 80 I = 1, K - 1
                  WORK(I) = WORK(I) + ABS(A(I,K))*XK
                  S = S + ABS(A(I,K))*ABS(X(I,J))
   80          CONTINUE
               WORK(K) = WORK(K) + ABS(A(K,K))*XK + S
  100       CONTINUE
         ELSE
            DO 140 K = 1, N
               S = ZERO
               XK = ABS(X(K,J))
               WORK(K) = WORK(K) + ABS(A(K,K))*XK
               DO 120 I = K + 1, N
                  WORK(I) = WORK(I) + ABS(A(I,K))*XK
                  S = S + ABS(A(I,K))*ABS(X(I,J))
  120          CONTINUE
               WORK(K) = WORK(K) + S
  140       CONTINUE
         END IF
         S = ZERO
         DO 160 I = 1, N
            IF (WORK(I).NE.ZERO) S = MAX(S,ABS(WORK(N+I))/WORK(I))
  160    CONTINUE
         BERR(J) = S
C
C        Test stopping criterion. Continue iterating if
C           1) The residual BERR(J) is larger than machine epsilon, and
C           2) BERR(J) decreased by at least a factor of 2 during the
C              last iteration, and
C           3) At most ITMAX iterations tried.
C
         IF (BERR(J).GT.EPS .AND. TWO*BERR(J)
     *       .LE.LSTRES .AND. COUNT.LE.ITMAX) THEN
C
C           Update solution and try again.
C
            CALL F07MEF(UPLO,N,1,AF,LDAF,IPIV,WORK(N+1),N,INFO)
            CALL DAXPY(N,ONE,WORK(N+1),1,X(1,J),1)
            LSTRES = BERR(J)
            COUNT = COUNT + 1
            GO TO 40
         END IF
C
C        Bound error from formula
C
C        norm(X - XTRUE) .le.
C        norm( abs(inv(A))*( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
C
C        where
C          norm(Z) is the magnitude of the largest component of Z
C          inv(A) is the inverse of A
C          abs(Z) is the componentwise absolute value of the matrix or
C             vector Z
C          NZ is the the maximum number of nonzeros in
C             any row of A, plus 1
C          EPS is machine epsilon
C
C        Use F04YCF to estimate the norm involving inv(A).
C
         DO 180 I = 1, N
            WORK(I) = ABS(WORK(N+I)) + NZ*EPS*WORK(I)
  180    CONTINUE
C
         KASE = 0
  200    CONTINUE
         IFAIL = 0
         CALL F04YCF(KASE,N,WORK(N+1),FERR(J),WORK(2*N+1),IWORK,IFAIL)
         IF (KASE.NE.0) THEN
            IF (KASE.EQ.1) THEN
C
C              Multiply by diag(R(1:N))*inv(A').
C
               CALL F07MEF(UPLO,N,1,AF,LDAF,IPIV,WORK(N+1),N,INFO)
               DO 220 I = 1, N
                  WORK(N+I) = WORK(I)*WORK(N+I)
  220          CONTINUE
            ELSE IF (KASE.EQ.2) THEN
C
C              Multiply by inv(A)*diag(R(1:N)).
C
               DO 240 I = 1, N
                  WORK(N+I) = WORK(I)*WORK(N+I)
  240          CONTINUE
               CALL F07MEF(UPLO,N,1,AF,LDAF,IPIV,WORK(N+1),N,INFO)
            END IF
            GO TO 200
         END IF
C
C        Normalize error.
C
         LSTRES = ZERO
         DO 260 I = 1, N
            LSTRES = MAX(LSTRES,ABS(X(I,J)))
  260    CONTINUE
         IF (LSTRES.NE.ZERO) FERR(J) = FERR(J)/LSTRES
C
  280 CONTINUE
C
      RETURN
C
C     End of F07MHF (DSYRFS)
C
      END
