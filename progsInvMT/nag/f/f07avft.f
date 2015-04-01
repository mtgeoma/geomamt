      SUBROUTINE F07AVF(TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX,
     *                  FERR,BERR,WORK,RWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZGERFS(TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,
     *                  LDX,FERR,BERR,WORK,RWORK,INFO)
C
C  Purpose
C  =======
C
C  ZGERFS improves the computed solution to a system of linear
C  equations and provides error bounds and backward error estimates for
C  the solutions.
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
C          of the matrices B and X.  NRHS >= 0.
C
C  A       (input) COMPLEX array, dimension (LDA,N)
C          The original n by n matrix A.
C
C  LDA     (input) INTEGER.
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  AF      (input) COMPLEX array, dimension (LDAF,N)
C          The factors L and U from the factorization A = P*L*U
C          as computed by F07ARF.
C
C  LDAF    (input) INTEGER
C          The leading dimension of the array AF.  LDAF >= max(1,N).
C
C  IPIV    (input) INTEGER array, dimension (N)
C          The pivot indices from F07ARF; for 1<=i<=N, row i of the
C          matrix was interchanged with row IPIV(i).
C
C  B       (input) COMPLEX array, dimension (LDB,NRHS)
C          The right hand side vectors for the system of linear
C          equations.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  X       (input/output) COMPLEX array, dimension (LDX,NRHS)
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
C  WORK    (workspace) COMPLEX array, dimension (2*N)
C
C  RWORK   (workspace) REAL array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C
C  Further Details
C  ===============
C
C  ITMAX is the maximum number of steps of iterative refinement.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      INTEGER           ITMAX
      PARAMETER         (ITMAX=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D+0)
      DOUBLE PRECISION  THREE
      PARAMETER         (THREE=3.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDAF, LDB, LDX, N, NRHS
      CHARACTER         TRANS
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), AF(LDAF,*), B(LDB,*), WORK(*),
     *                  X(LDX,*)
      DOUBLE PRECISION  BERR(*), FERR(*), RWORK(*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      COMPLEX*16        ZDUM
      DOUBLE PRECISION  EPS, LSTRES, S, XK
      INTEGER           COUNT, I, IFAIL, J, K, KASE, NZ
      LOGICAL           NOTRAN
      CHARACTER         TRANSN, TRANST
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          ZAXPY, ZCOPY, ZGEMV, F04ZCF, F06AAZ, F07ASF
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
      ELSE IF (LDAF.LT.MAX(1,N)) THEN
         INFO = -7
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -10
      ELSE IF (LDX.LT.MAX(1,N)) THEN
         INFO = -12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07AVF/ZGERFS',-INFO)
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
      IF (NOTRAN) THEN
         TRANSN = 'N'
         TRANST = 'C'
      ELSE
         TRANSN = 'C'
         TRANST = 'N'
      END IF
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
C        Compute residual R = B - op(A) * X,
C        where op(A) = A, A**T, or A**H, depending on TRANS.
C
         CALL ZCOPY(N,B(1,J),1,WORK,1)
         CALL ZGEMV(TRANS,N,N,-ONE,A,LDA,X(1,J),1,ONE,WORK,1)
C
C        Compute componentwise relative backward error from formula
C
C        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
C
C        where 0/0 is treated as 0, and abs(Z) is the componentwise
C        absolute value of the matrix or vector Z.
C
         DO 60 I = 1, N
            RWORK(I) = CABS1(B(I,J))
   60    CONTINUE
C
C        Compute abs(op(A))*abs(X) + abs(B).
C
         IF (NOTRAN) THEN
            DO 100 K = 1, N
               XK = CABS1(X(K,J))
               DO 80 I = 1, N
                  RWORK(I) = RWORK(I) + CABS1(A(I,K))*XK
   80          CONTINUE
  100       CONTINUE
         ELSE
            DO 140 K = 1, N
               S = ZERO
               DO 120 I = 1, N
                  S = S + CABS1(A(I,K))*CABS1(X(I,J))
  120          CONTINUE
               RWORK(K) = RWORK(K) + S
  140       CONTINUE
         END IF
         S = ZERO
         DO 160 I = 1, N
            IF (RWORK(I).NE.ZERO) S = MAX(S,CABS1(WORK(I))/RWORK(I))
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
            CALL F07ASF(TRANS,N,1,AF,LDAF,IPIV,WORK,N,INFO)
            CALL ZAXPY(N,ONE,WORK,1,X(1,J),1)
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
C        Note A is replaced by A**H if TRANS = 'T' or 'C'.
C
C        Use F04ZCF to estimate the infinity-norm of the matrix
C           inv(op(A)) * diag(W),
C        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
C
         DO 180 I = 1, N
            RWORK(I) = CABS1(WORK(I)) + NZ*EPS*RWORK(I)
  180    CONTINUE
C
         KASE = 0
  200    CONTINUE
         IFAIL = 0
         CALL F04ZCF(KASE,N,WORK,FERR(J),WORK(N+1),IFAIL)
         IF (KASE.NE.0) THEN
            IF (KASE.EQ.1) THEN
C
C              Multiply by diag(W)*inv(op(A)**H).
C
               CALL F07ASF(TRANST,N,1,AF,LDAF,IPIV,WORK,N,INFO)
               DO 220 I = 1, N
                  WORK(I) = RWORK(I)*WORK(I)
  220          CONTINUE
            ELSE
C
C              Multiply by inv(op(A))*diag(W).
C
               DO 240 I = 1, N
                  WORK(I) = RWORK(I)*WORK(I)
  240          CONTINUE
               CALL F07ASF(TRANSN,N,1,AF,LDAF,IPIV,WORK,N,INFO)
            END IF
            GO TO 200
         END IF
C
C        Normalize error.
C
         LSTRES = ZERO
         DO 260 I = 1, N
            LSTRES = MAX(LSTRES,CABS1(X(I,J)))
  260    CONTINUE
         IF (LSTRES.NE.ZERO) FERR(J) = FERR(J)/LSTRES
C
  280 CONTINUE
C
      RETURN
C
C     End of F07AVF (ZGERFS)
C
      END
