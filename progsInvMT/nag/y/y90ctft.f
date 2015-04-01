      SUBROUTINE Y90CTF(MATRIX,JOBS,N,NJORD,JORD,EIG,A,LDA,COND,X,LDX,Y,
     *                  LDY,Q,LDQ,VEC,WK1,LDWK1,SEED)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C-----------------------------------------------------------------------
C
C         ==================================================
C         *  Y90CTF :  Unsymmetric Eigenproblem Generator  *
C         ==================================================
C
C
C     DESCRIPTION
C     ===========
C
C     Y90RJF generates a matrix of known eigenvalues and eigenvectors.
C     It accepts in input the structure of the Jordan matrix, such,
C     that is, that the matrix A is given by:
C
C              A  =  X * J * X**(-1)
C
C     where X is non-singular.  The condition number of X can also be
C     determined in input.
C
C     Two points need to be made:
C
C     1. Eigenvalue or cluster condition numbers depend on the condition
C        number of X.
C     2. Eigenvectors and invariant subspace condition numbers depend
C        the condition number of X, the spacing between eigenvalues and
C        the structure of the Jordan matrix.
C
C     In output the following items can be returned:
C
C     1. A full square matrix A or a real Schur matrix S.
C     2. The matrices X and X**(-1) of, respectively, the right- and
C        left-eigenvectorsof A or S.
C     3. The matrix of the Schur vectors Q.
C
C
C     DEPENDENCIES
C     ============
C
C     Y90CTF depend on the prior successful testing of the following
C     LAPACK routines:
C
C     1. DGEQRF (QR factorization)
C     2. ZUNGQR (Computation of the orthogonal matrix Q)
C     3. ZUNMQR (Multiplication of a matrix by Q)
C
C
C     ARGUMENTS
C     =========
C
C     MATRIX :  Character*1, input.
C               Defines the type of matrix required:
C               'F' or 'f'  ==>  Full square matrix.
C               'S' or 's'  ==>  Schur matrix.
C
C     JOBS   :  Character*1, input.
C               Defines the Schur vectors computations requests:
C               'N' or 'n'  ==>  No Schur vectors required.
C               'V' or 'v'  ==>  Schur vectors.
C
C     N      :  Integer, input.
C               Order of the matrices to be generated.
C
C     NJORD  :  Integer, input.
C               Number of Jordan blocks.
C
C     JORD   :  Integer array of DIMENSION (*), input.
C               The absolute value of JORD(i) denotes the size of the
C               i-th Jordan block.  For complex eigenpairs, JORD(i)
C               is negative and only one of the pair needs to be
C               specified, and only one Jordan block.
C               It must be:
C                  SUM(j=1,...,NJORD) ABS(JORD(i))*c(i) = N
C               where c(i) is 1 for real eigenvalues, 2 for complex
C               eigenvalues.
C
C     EIG    :  Complex array of DIMENSION (*), input/output.
C               On input:  EIG(i) contains the eigenvalue for the
C                          i-th Jordan block (NJORD values in all).
C               On output: EIG(i) contains the i-th eigenvalue (N
C                          eigenvalues in all).
C               Naturally, different Jordan blocks may have identical
C               eigenvalues.
C
C     A      :  Complex array of DIMENSION (LDA,*), output.
C               Contains the square or Schur matrix to be generated.
C               It must have at least N columns.
C
C     LDA    :  Integer, input.
C               Leading dimension of the array A.  It must be LDA >= N.
C
C     COND   :  Complex, input.
C               Condition number for the matrix X to be generated.
C               Singular values in geometric progression are generated,
C               such that sigma(max)/sigma(min) = COND.
C
C     X      :  Complex array of DIMENSION (LDX,*), output.
C               Contains the matrix X. It must have at least N columns.
C
C     LDX    :  Integer, input.
C               Leading dimension of the array X.  It must be LDX >= N.
C
C     Y     :  Complex array of DIMENSION (LDY,*), output.
C               Contains the matrix X**(-1). It must have at least
C               N columns.
C
C     LDY   :  Integer, input.
C               Leading dimension of the array Y.
C               It must be LDY >= N
C
C     Q      :  Complex array of DIMENSION (LDQ,*), output.
C               Contains the matrix Q of the Schur vectors.  It must
C               have at least N columns.
C               NOTE: Q is also used as workspace.
C
C     LDQ    :  Integer, input.
C               Leading dimension of the array Q.  It must be LDQ >= N.
C
C     VEC    :  Complex array of DIMENSION (*), workspace.
C
C     WK1    :  Complex array of DIMENSION (LDWK1,*), workspace.
C               It must have at least N columns.
C
C     LDWK1  :  Integer, input.
C               Leading dimension of WK1.  It must be LDWK1 >= N.
C
C     SEED   :  Integer array of DIMENSION (4), input/output.
C               Seeds used by the random number generator.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND
      INTEGER           LDA, LDQ, LDWK1, LDX, LDY, N, NJORD
      CHARACTER*1       JOBS, MATRIX
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,N), EIG(N), Q(LDQ,N), VEC(N),
     *                  WK1(LDWK1,N), X(LDX,N), Y(LDY,N)
      INTEGER           JORD(N), SEED(4)
C     .. Local Scalars ..
      COMPLEX*16        DETMAN, SCALE
      INTEGER           DETEXP, DIST, DTYPE, I, INFO, J, K, K1
C     .. Local Arrays ..
      COMPLEX*16        DIAG(2)
C     .. External Subroutines ..
      EXTERNAL          F06GCF, F06GDF, F06HBF, F06THF, F06ZAF, Y90CDF,
     *                  Y90CGF, Y90DHF, ZGEQRF, ZUNGQR, ZUNMQR
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize
C
C-----------------------------------------------------------------------
      DIST = 1
      DTYPE = 2
      SCALE = CONE
C-----------------------------------------------------------------------
C
C     1. Define the eigenvalues and the Jordan structure
C
C-----------------------------------------------------------------------
      K1 = N + 1
      DO 20 I = NJORD, 1, -1
         K1 = K1 - JORD(I)
         CALL F06HBF(JORD(I),EIG(I),EIG(K1),1)
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     2. Generate the singular values of the matrix X and two random
C        orthogonal matrices
C
C-----------------------------------------------------------------------
C
C     2.1 Generate the singular values
C
      CALL Y90CGF('N',DTYPE,N,VEC,DIAG,COND,SCALE,DETMAN,DETEXP,DIST,
     *            SEED)
C
C     2.2 Generate the random orthogonal matrices
C
      CALL Y90CDF('Left','I',N,N,WK1,LDWK1,X,Y,SEED)
      CALL Y90CDF('Right','I',N,N,A,LDA,X,Y,SEED)
C-----------------------------------------------------------------------
C
C     3. Generate the right eigenvectors
C
C-----------------------------------------------------------------------
      CALL Y90DHF('N',N,N,WK1,LDWK1,Q,LDQ)
      DO 40 I = 1, N
         CALL F06GDF(N,VEC(I),Q(1,I),1)
   40 CONTINUE
      CALL F06ZAF('N','C',N,N,N,CONE,Q,LDQ,A,LDA,CZERO,X,LDX)
C-----------------------------------------------------------------------
C
C     4. Generate the left eigenvectors
C
C-----------------------------------------------------------------------
      DO 60 I = 1, N
         CALL F06GDF(N,CONE/VEC(I),WK1(1,I),1)
   60 CONTINUE
      CALL F06ZAF('N','C',N,N,N,CONE,WK1,LDWK1,A,LDA,CZERO,Y,LDY)
C-----------------------------------------------------------------------
C
C     5. In order to generate a Schur matrix, the matrix of the right
C        eigenvectors is reduced to upper triangular form, and the same
C        transformation is then applied to the matrix of the left
C        eigenvectors.
C
C-----------------------------------------------------------------------
      IF ((MATRIX.EQ.'S') .OR. (MATRIX.EQ.'s')) THEN
C
C     5.1 Carry ou the QR factorization of X to obtain the right
C         eigenvectors of the Schur matrix.
C
         CALL ZGEQRF(N,N,X,LDX,VEC,WK1,N*LDWK1,INFO)
C
C     5.2 Multiply Y by Q (the orthogonal matrix of the QR
C         factorization) to obtain the left eigenvectors of the Schur
C         matrix.
C
         CALL ZUNMQR('L','C',N,N,N,X,LDX,VEC,Y,LDY,WK1,N*LDWK1,INFO)
         CALL F06THF('L',N-1,N-1,CZERO,CZERO,X(2,1),LDX)
         CALL F06THF('U',N-1,N-1,CZERO,CZERO,Y(1,2),LDY)
C
C     5.3 Compute the matrix of the Schur vectors, if so required
C
         IF ((JOBS.EQ.'V' .OR. JOBS.EQ.'v')) THEN
            CALL Y90DHF('N',N,N,X,LDX,Q,LDQ)
            CALL ZUNGQR(N,N,N,Q,LDQ,VEC,WK1,N*LDWK1,INFO)
         END IF
C
      END IF
C-----------------------------------------------------------------------
C
C     6. Compute the matrix
C
C-----------------------------------------------------------------------
      CALL Y90DHF('N',N,N,X,LDX,WK1,LDWK1)
      DO 80 I = 1, N
         CALL F06GDF(N,EIG(I),WK1(1,I),1)
   80 CONTINUE
C
      K = 1
      DO 120 I = 1, NJORD
         DO 100 J = K, K + JORD(I) - 2
            CALL F06GCF(N,CONE,X(1,J),1,WK1(1,J+1),1)
  100    CONTINUE
         K = K + JORD(I)
  120 CONTINUE
C
      CALL F06ZAF('N','C',N,N,N,CONE,WK1,LDWK1,Y,LDY,CZERO,A,LDA)
C-----------------------------------------------------------------------
C
C     End of subroutine Y90CTF
C
C-----------------------------------------------------------------------
      RETURN
      END
