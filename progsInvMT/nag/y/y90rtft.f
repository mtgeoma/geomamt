      SUBROUTINE Y90RTF(MATRIX,JOBS,N,NJORD,JORD,EIGR,EIGI,A,LDA,COND,X,
     *                  LDX,YH,LDYH,Q,LDQ,VEC,IVEC,WK1,LDWK1,SEED)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C-----------------------------------------------------------------------
C
C         ==================================================
C         *  Y90RTF :  Unsymmetric Eigenproblem Generator  *
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
C     Y90RTF depend on the prior successful testing of the following
C     LAPACK routines:
C
C     1. DGEQRF (QR factorization)
C     2. DORGQR (Computation of the orthogonal matrix Q)
C     3. DORMQR (Multiplication of a matrix by Q)
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
C     EIGR   :  Real array of DIMENSION (*), input/output.
C               On input:  EIGR(i) contains the real part of the
C                          eigenvalue for the i-th Jordan block
C                          (NJORD values in all).
C               On output: EIGR(i) contains the real part of the
C                          i-th eigenvalue (N eigenvalues in all).
C               Notice that only one Jordan block for each conjugate
C               pair, hence only one value of EIGR and EIGI (see below)
C               need to be specified.  Naturally, different Jordan
C               blocks may have identical eigenvalues.
C
C     EIGI   :  Real array of DIMENSION (*), input/output.
C               Same as EIGR for the imaginary part of the eigenvalues.
C
C     A      :  Real array of DIMENSION (LDA,*), output.
C               Contains the square or Schur matrix to be generated.
C               It must have at least N columns.
C
C     LDA    :  Integer, input.
C               Leading dimension of the array A.  It must be LDA >= N.
C
C     COND   :  Real, input.
C               Condition number for the matrix X to be generated.
C               Singular values in geometric progression are generated,
C               such that sigma(max)/sigma(min) = COND.
C
C     X      :  Real array of DIMENSION (LDX,*), output.
C               Contains the matrix X. It must have at least N columns.
C
C     LDX    :  Integer, input.
C               Leading dimension of the array X.  It must be LDX >= N.
C
C     YH     :  Real array of DIMENSION (LDYH,*), output.
C               Contains the matrix X**(-1). It must have at least
C               N columns.
C
C     LDYH   :  Integer, input.
C               Leading dimension of the array YH.
C               It must be LDYH >= N
C
C     Q      :  Real array of DIMENSION (LDQ,*), output.
C               Contains the matrix Q of the Schur vectors.  It must
C               have at least N columns.
C               NOTE: Q is also used as workspace.
C
C     LDQ    :  Integer, input.
C               Leading dimension of the array Q.  It must be LDQ >= N.
C
C     VEC    :  Real array of DIMENSION (*), workspace.
C
C     IVEC   :  Integer array of DIMENSION (*), workspace.
C
C     WK1    :  Real array of DIMENSION (LDWK1,*), workspace.
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
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND
      INTEGER           LDA, LDQ, LDWK1, LDX, LDYH, N, NJORD
      CHARACTER*1       JOBS, MATRIX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), EIGI(*), EIGR(*), Q(LDQ,*), VEC(*),
     *                  WK1(LDWK1,*), X(LDX,*), YH(LDYH,*)
      INTEGER           IVEC(*), JORD(*), SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  CN, DETMAN, P, SCALE, SIGMA, SN, TAU, TEMP
      INTEGER           DETEXP, DIST, DTYPE, I, INFO, J, K, K1, NCJ, NJ
      LOGICAL           PAIR
C     .. Local Arrays ..
      DOUBLE PRECISION  DIAG(2)
C     .. External Subroutines ..
      EXTERNAL          DGEQRF, DORGQR, DORMQR, F06EDF, F06EGF, F06EPF,
     *                  F06FBF, F06QHF, F06YAF, Y90RDF, Y90RGF, Y90RTX,
     *                  Y90SHF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize
C
C-----------------------------------------------------------------------
      DIST = 1
      DTYPE = 2
      SCALE = ONE
C-----------------------------------------------------------------------
C
C     1. Define the eigenvalues and the Jordan structure
C
C-----------------------------------------------------------------------
      K1 = N + 1
      DO 20 I = NJORD, 1, -1
         IF (JORD(I).GT.0) THEN
            K1 = K1 - JORD(I)
            CALL F06FBF(JORD(I),EIGR(I),EIGR(K1),1)
            CALL F06FBF(JORD(I),ZERO,EIGI(K1),1)
         ELSE
            K1 = K1 + 2*JORD(I)
            CALL F06FBF(-2*JORD(I),EIGR(I),EIGR(K1),1)
            CALL F06FBF(-JORD(I),EIGI(I),EIGI(K1),2)
            CALL F06FBF(-JORD(I),-EIGI(I),EIGI(K1+1),2)
         END IF
   20 CONTINUE
C
      NCJ = 0
      PAIR = .FALSE.
      DO 40 I = 1, N
         IF (PAIR) THEN
            PAIR = .FALSE.
         ELSE
            IF (EIGI(I).NE.ZERO) THEN
               PAIR = .TRUE.
               NCJ = NCJ + 1
               IVEC(NCJ) = I
            END IF
         END IF
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     2. Generate the singular values of the matrix X and two random
C        orthogonal matrices
C
C-----------------------------------------------------------------------
C
C     2.1 Generate the singular values
C
      NJ = N - NCJ
      CALL Y90RGF('N',DTYPE,NJ,VEC,DIAG,COND,SCALE,DETMAN,DETEXP,DIST,
     *            SEED)
C
      DO 80 I = 1, NCJ
         K = IVEC(I)
         DO 60 J = NJ, K + 1, -1
            VEC(J+1) = VEC(J)
   60    CONTINUE
         NJ = NJ + 1
         VEC(K+1) = VEC(K)
   80 CONTINUE
C
C     2.2 Generate the random orthogonal matrices
C
      CALL Y90RDF('Left','I',N,N,WK1,LDWK1,X,YH,SEED)
      CALL Y90RDF('Right','I',N,N,A,LDA,X,YH,SEED)
C-----------------------------------------------------------------------
C
C     3. Generate the right eigenvectors
C
C-----------------------------------------------------------------------
      CALL Y90SHF('N',N,N,WK1,LDWK1,Q,LDQ)
      DO 100 I = 1, N
         CALL F06EDF(N,VEC(I),Q(1,I),1)
  100 CONTINUE
      CALL F06YAF('N','T',N,N,N,ONE,Q,LDQ,A,LDA,ZERO,X,LDX)
C-----------------------------------------------------------------------
C
C     4. Generate the left eigenvectors
C
C-----------------------------------------------------------------------
      DO 120 I = 1, N
         CALL F06EDF(N,ONE/VEC(I),WK1(1,I),1)
  120 CONTINUE
      CALL F06YAF('N','T',N,N,N,ONE,WK1,LDWK1,A,LDA,ZERO,YH,LDYH)
C-----------------------------------------------------------------------
C
C     5. Generate a Schur matrix, the matrices of its right and left
C        eigenvectors and the matrix of the Schur vectors if required
C
C-----------------------------------------------------------------------
      IF ((MATRIX.EQ.'S') .OR. (MATRIX.EQ.'s')) THEN
C
C     5.1 Carry ou the QR factorization of X to obtain the right
C         eigenvectors of the Schur matrix.
C
         CALL DGEQRF(N,N,X,LDX,VEC,WK1,N*LDWK1,INFO)
C
C     5.2 Multiply YH by Q (the orthogonal matrix of the QR
C         factorization) to obtain the left eigenvectors of the Schur
C         matrix.
C
         CALL DORMQR('L','T',N,N,N,X,LDX,VEC,YH,LDYH,WK1,N*LDWK1,INFO)
         CALL F06QHF('L',N-1,N-1,ZERO,ZERO,X(2,1),LDX)
         CALL F06QHF('U',N-1,N-1,ZERO,ZERO,YH(1,2),LDYH)
C
C     5.3 Compute the matrix of the Schur vectors, if so required
C
         IF ((JOBS.EQ.'V' .OR. JOBS.EQ.'v')) THEN
            CALL Y90SHF('N',N,N,X,LDX,Q,LDQ)
            CALL DORGQR(N,N,N,Q,LDQ,VEC,WK1,N*LDWK1,INFO)
         END IF
C
C     Generate the Schur matrix (NOTE: This is NOT in canonical form!)
C
         CALL Y90RTX(N,NJORD,JORD,EIGR,EIGI,A,LDA,X,LDX,YH,LDYH,WK1,
     *               LDWK1)
C
C     Compute the transformations that reduce the Schur matrix to
C     canonical form
C
         DO 140 I = 1, NCJ
            K = IVEC(I)
            TEMP = A(K,K) - A(K+1,K+1)
            P = HALF*TEMP
            SIGMA = A(K,K+1) + A(K+1,K)
            TAU = SQRT(SIGMA*SIGMA+TEMP*TEMP)
            IF (TAU.NE.ZERO) THEN
               CN = SQRT(HALF*(ONE+ABS(SIGMA)/TAU))
               SN = -(P/(TAU*CN))*SIGN(ONE,SIGMA)
               CALL F06EPF(N,X(K,1),LDX,X(K+1,1),LDX,CN,SN)
               CALL F06EPF(N,YH(K,1),LDYH,YH(K+1,1),LDYH,CN,SN)
               CALL F06EPF(N,A(K,1),LDA,A(K+1,1),LDA,CN,SN)
               CALL F06EPF(N,A(1,K),1,A(1,K+1),1,CN,SN)
            END IF
            IF (A(K+1,K).LT.ZERO) THEN
               CALL F06EGF(N,A(1,K),1,A(1,K+1),1)
               CALL F06EGF(N,A(K,1),LDA,A(K+1,1),LDA)
               CALL F06EGF(N,X(K,1),LDX,X(K+1,1),LDX)
               CALL F06EGF(N,YH(K,1),LDYH,YH(K+1,1),LDYH)
            END IF
  140    CONTINUE
C-----------------------------------------------------------------------
C
C     6. Generate a full square matrix
C
C-----------------------------------------------------------------------
      ELSE
C
         CALL Y90RTX(N,NJORD,JORD,EIGR,EIGI,A,LDA,X,LDX,YH,LDYH,WK1,
     *               LDWK1)
C
      END IF
C-----------------------------------------------------------------------
C
C     7. Normalize to unity the complex eigenvectors.
C
C-----------------------------------------------------------------------
      SCALE = SQRT(HALF)
      DO 160 I = 1, NCJ
         K = IVEC(I)
         CALL F06EDF(N,SCALE,X(1,K),1)
         CALL F06EDF(N,SCALE,X(1,K+1),1)
         CALL F06EDF(N,SCALE,YH(1,K),1)
         CALL F06EDF(N,SCALE,YH(1,K+1),1)
  160 CONTINUE
C-----------------------------------------------------------------------
C
C     End of subroutine Y90RTF
C
C-----------------------------------------------------------------------
      RETURN
      END
