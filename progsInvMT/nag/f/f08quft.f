      SUBROUTINE F08QUF(JOB,COMPQ,SELECT,N,T,LDT,Q,LDQ,W,M,S,SEP,WORK,
     *                  LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZTRSEN(JOB,COMPQ,SELECT,N,T,LDT,Q,LDQ,W,M,S,SEP,
     *                  WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZTRSEN reorders the Schur factorization of a complex matrix
C  A = Q*T*Q', so that a selected cluster of eigenvalues appears in the
C  leading positions on the diagonal of the upper triangular matrix T,
C  and the leading columns of Q form an orthonormal basis of the
C  corresponding right invariant subspace.
C
C  Optionally the routine computes the reciprocal condition numbers of
C  the cluster of eigenvalues and/or the invariant subspace.
C
C  Arguments
C  =========
C
C  JOB     (input) CHARACTER*1
C          Specifies whether condition numbers are required for the
C          cluster of eigenvalues (S) or the invariant subspace (SEP):
C          = 'N': none;
C          = 'E': for eigenvalues only (S);
C          = 'V': for invariant subspace only (SEP);
C          = 'B': for both eigenvalues and invariant subspace (S and
C                 SEP).
C
C  COMPQ   (input) CHARACTER*1
C          = 'V': update the matrix Q of Schur vectors;
C          = 'N': do not update Q.
C
C  SELECT  (input) LOGICAL array, dimension (N)
C          SELECT specifies the eigenvalues in the selected cluster. To
C          select the j-th eigenvalue, SELECT(j) must be set to .TRUE..
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input/output) COMPLEX*16 array, dimension(LDT,N)
C          On entry, the upper triangular matrix T.
C          On exit, T is overwritten by the reordered matrix T, with the
C          selected eigenvalues in the leading positions on the
C          diagonal.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)
C          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
C          On exit, if COMPQ = 'V', Q has been postmultiplied by the
C          unitary transformation matrix which reorders T; the leading M
C          columns of Q form an orthonormal basis for the specified
C          invariant subspace.
C          If COMPQ = 'N', Q is not referenced.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q.
C          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
C
C  W       (output) COMPLEX*16
C          The reordered eigenvalues of T, in the same order as they
C          appear on the diagonal of the output matrix T.
C
C  M       (output) INTEGER
C          The dimension of the specified invariant subspace.
C
C  S       (output) DOUBLE PRECISION
C          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
C          condition number for the selected cluster of eigenvalues.
C          S cannot underestimate the true reciprocal condition number
C          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
C          If JOB = 'N' or 'V', S is not referenced.
C
C  SEP     (output) DOUBLE PRECISION
C          If JOB = 'V' or 'B', SEP is the estimated reciprocal
C          condition number of the specified invariant subspace. If
C          M = 0 or N, SEP = norm(T).
C          If JOB = 'N' or 'E', SEP is not referenced.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          If JOB = 'N', WORK is not referenced.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.
C          If JOB = 'N', LWORK >= 1;
C          if JOB = 'E', LWORK = M*(N-M);
C          if JOB = 'V' or 'B', LWORK >= 2*M*(N-M).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  Further Details
C  ===============
C
C  ZTRSEN first collects the selected eigenvalues by computing a unitary
C  transformation Z to move them to the top left corner of T. In other
C  words, the selected eigenvalues are the eigenvalues of T11 in:
C
C                Z'*T*Z = ( T11 T12 ) n1
C                         (  0  T22 ) n2
C                            n1  n2
C
C  where N = n1+n2 and Z' means the conjugate transpose of Z. The first
C  n1 columns of Z span the specified invariant subspace of T.
C
C  If T has been obtained from the Schur factorization of a matrix
C  A = Q*T*Q', then the reordered Schur factorization of A is given by
C  A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span the
C  corresponding invariant subspace of A.
C
C  The reciprocal condition number of the average of the eigenvalues of
C  T11 may be returned in S. S lies between 0 (very badly conditioned)
C  and 1 (very well conditioned). It is computed as follows. First we
C  compute R so that
C
C                         P = ( I  R ) n1
C                             ( 0  0 ) n2
C                               n1 n2
C
C  is the projector on the invariant subspace associated with T11.
C  R is the solution of the Sylvester equation:
C
C                        T11*R - R*T22 = T12.
C
C  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
C  the two-norm of M. Then S is computed as the lower bound
C
C                      (1 + F-norm(R)**2)**(-1/2)
C
C  on the reciprocal of 2-norm(P), the true reciprocal condition number.
C  S cannot underestimate 1 / 2-norm(P) by more than a factor of
C  sqrt(N).
C
C  An approximate error bound for the computed average of the
C  eigenvalues of T11 is
C
C                         EPS * norm(T) / S
C
C  where EPS is the machine precision.
C
C  The reciprocal condition number of the right invariant subspace
C  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
C  SEP is defined as the separation of T11 and T22:
C
C                     sep( T11, T22 ) = sigma-min( C )
C
C  where sigma-min(C) is the smallest singular value of the
C  n1*n2-by-n1*n2 matrix
C
C     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
C
C  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
C  product. We estimate sigma-min(C) by the reciprocal of an estimate of
C  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
C  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
C
C  When SEP is small, small changes in T can cause large changes in
C  the invariant subspace. An approximate bound on the maximum angular
C  error in the computed right invariant subspace is
C
C                      EPS * norm(T) / SEP
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S, SEP
      INTEGER           INFO, LDQ, LDT, LWORK, M, N
      CHARACTER         COMPQ, JOB
C     .. Array Arguments ..
      COMPLEX*16        Q(LDQ,*), T(LDT,*), W(*), WORK(LWORK)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EST, RNORM, SCALE
      INTEGER           IERR, IFAIL, K, KASE, KS, N1, N2, NN
      LOGICAL           WANTBH, WANTQ, WANTS, WANTSP
C     .. Local Arrays ..
      DOUBLE PRECISION  RWORK(1)
C     .. External Functions ..
      DOUBLE PRECISION  F06UAF
      EXTERNAL          F06UAF
C     .. External Subroutines ..
      EXTERNAL          F04ZCF, F06AAZ, F06TFF, ZTREXC, ZTRSYL
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, SQRT
C     .. Executable Statements ..
C
C     Decode and test the input parameters.
C
      WANTBH = (JOB.EQ.'B' .OR. JOB.EQ.'b')
      WANTS = (JOB.EQ.'E' .OR. JOB.EQ.'e') .OR. WANTBH
      WANTSP = (JOB.EQ.'V' .OR. JOB.EQ.'v') .OR. WANTBH
      WANTQ = (COMPQ.EQ.'V' .OR. COMPQ.EQ.'v')
C
C     Set M to the number of selected eigenvalues.
C
      M = 0
      DO 20 K = 1, N
         IF (SELECT(K)) M = M + 1
   20 CONTINUE
C
      N1 = M
      N2 = N - M
      NN = N1*N2
C
      INFO = 0
      IF ( .NOT. (JOB.EQ.'N' .OR. JOB.EQ.'n')
     *    .AND. .NOT. WANTS .AND. .NOT. WANTSP) THEN
         INFO = -1
      ELSE IF ( .NOT. (COMPQ.EQ.'N' .OR. COMPQ.EQ.'n')
     *         .AND. .NOT. WANTQ) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (LDT.LT.MAX(1,N)) THEN
         INFO = -6
      ELSE IF (LDQ.LT.1 .OR. (WANTQ .AND. LDQ.LT.N)) THEN
         INFO = -8
      ELSE IF (LWORK.LT.1 .OR. ((WANTS .AND. .NOT. WANTSP)
     *         .AND. LWORK.LT.NN) .OR. (WANTSP .AND. LWORK.LT.2*NN))
     *         THEN
         INFO = -14
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QUF/ZTRSEN',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.N .OR. M.EQ.0) THEN
         IF (WANTS) S = ONE
         IF (WANTSP) SEP = F06UAF('1',N,N,T,LDT,RWORK)
         GO TO 80
      END IF
C
C     Collect the selected eigenvalues at the top left corner of T.
C
      KS = 0
      DO 40 K = 1, N
         IF (SELECT(K)) THEN
            KS = KS + 1
C
C           Swap the K-th eigenvalue to position KS.
C
            IF (K.NE.KS) CALL ZTREXC(COMPQ,N,T,LDT,Q,LDQ,K,KS,IERR)
         END IF
   40 CONTINUE
C
      IF (WANTS) THEN
C
C        Solve the Sylvester equation for R:
C
C           T11*R - R*T22 = scale*T12
C
         CALL F06TFF('General',N1,N2,T(1,N1+1),LDT,WORK,N1)
         CALL ZTRSYL('N','N',-1,N1,N2,T,LDT,T(N1+1,N1+1),LDT,WORK,N1,
     *               SCALE,IERR)
C
C        Estimate the reciprocal of the condition number of the cluster
C        of eigenvalues.
C
         RNORM = F06UAF('F',N1,N2,WORK,N1,RWORK)
         IF (RNORM.EQ.ZERO) THEN
            S = ONE
         ELSE
            S = SCALE/(SQRT(SCALE*SCALE/RNORM+RNORM)*SQRT(RNORM))
         END IF
      END IF
C
      IF (WANTSP) THEN
C
C        Estimate sep(T11,T22).
C
         EST = ZERO
         KASE = 0
   60    CONTINUE
         IFAIL = 0
         CALL F04ZCF(KASE,NN,WORK,EST,WORK(NN+1),IFAIL)
         IF (KASE.NE.0) THEN
            IF (KASE.EQ.1) THEN
C
C              Solve T11*R - R*T22 = scale*X.
C
               CALL ZTRSYL('N','N',-1,N1,N2,T,LDT,T(N1+1,N1+1),LDT,WORK,
     *                     N1,SCALE,IERR)
            ELSE
C
C              Solve T11'*R - R*T22' = scale*X.
C
               CALL ZTRSYL('C','C',-1,N1,N2,T,LDT,T(N1+1,N1+1),LDT,WORK,
     *                     N1,SCALE,IERR)
            END IF
            GO TO 60
         END IF
C
         SEP = SCALE/EST
      END IF
C
   80 CONTINUE
C
C     Copy reordered eigenvalues to W.
C
      DO 100 K = 1, N
         W(K) = T(K,K)
  100 CONTINUE
      RETURN
C
C     End of F08QUF (ZTRSEN)
C
      END
