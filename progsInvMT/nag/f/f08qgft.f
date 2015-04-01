      SUBROUTINE F08QGF(JOB,COMPQ,SELECT,N,T,LDT,Q,LDQ,WR,WI,M,S,SEP,
     *                  WORK,LWORK,IWORK,LIWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DTRSEN(JOB,COMPQ,SELECT,N,T,LDT,Q,LDQ,WR,WI,M,S,
     *                  SEP,WORK,LWORK,IWORK,LIWORK,INFO)
C
C  Purpose
C  =======
C
C  DTRSEN reorders the real Schur factorization of a real matrix
C  A = Q*T*Q', so that a selected cluster of eigenvalues appears in the
C  leading diagonal blocks of the upper quasi-triangular matrix T,
C  and the leading columns of Q form an orthonormal basis of the
C  corresponding right invariant subspace.
C
C  Optionally the routine computes the reciprocal condition numbers of
C  the cluster of eigenvalues and/or the invariant subspace.
C
C  T must be in Schur canonical form, that is, block upper triangular
C  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
C  has its diagonal elemnts equal and its off-diagonal elements of
C  opposite sign.
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
C          select a real eigenvalue w(j), SELECT(j) must be set to
C          .TRUE.. To select a complex conjugate pair of eigenvalues
C          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
C          either SELECT(j) or SELECT(j+1) or both must be set to
C          .TRUE.; a complex conjugate pair of eigenvalues must be
C          either both included in the cluster or both excluded.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input/output) DOUBLE PRECISION array, dimension(LDT,N)
C          On entry, the upper quasi-triangular matrix T, in Schur
C          canonical form.
C          On exit, T is overwritten by the reordered matrix T, again in
C          Schur canonical form, with the selected eigenvalues in the
C          leading diagonal blocks.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
C          On exit, if COMPQ = 'V', Q has been postmultiplied by the
C          orthogonal transformation matrix which reorders T; the
C          leading M columns of Q form an orthonormal basis for the
C          specified invariant subspace.
C          If COMPQ = 'N', Q is not referenced.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q.
C          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
C
C  WR      (output) DOUBLE PRECISION array, dimension (N)
C  WI      (output) DOUBLE PRECISION array, dimension (N)
C          The real and imaginary parts, respectively, of the reordered
C          eigenvalues of T. The eigenvalues are stored in the same
C          order as on the diagonal of T, with WR(i) = T(i,i) and, if
C          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
C          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
C          sufficiently ill-conditioned, then its value may differ
C          significantly from its value before reordering.
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
C  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.
C          If JOB = 'N', LWORK >= max(1,N);
C          if JOB = 'E', LWORK >= M*(N-M);
C          if JOB = 'V' or 'B', LWORK >= 2*M*(N-M).
C
C  IWORK   (workspace) INTEGER
C          IF JOB = 'N' or 'E', IWORK is not referenced.
C
C  LIWORK  (input) INTEGER
C          The dimension of the array IWORK.
C          If JOB = 'N' or 'E', LIWORK >= 1;
C          if JOB = 'V' or 'B', LIWORK >= M*(N-M).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C          = 1: reordering of T failed because some eigenvalues are too
C               close to separate (the problem is very ill-conditioned);
C               T may have been partially reordered, and WR and WI
C               contain the eigenvalues in the same order as in T; S and
C               SEP (if requested) are set to zero.
C
C  Further Details
C  ===============
C
C  DTRSEN first collects the selected eigenvalues by computing an
C  orthogonal transformation Z to move them to the top left corner of T.
C  In other words, the selected eigenvalues are the eigenvalues of T11
C  in:
C
C                Z'*T*Z = ( T11 T12 ) n1
C                         (  0  T22 ) n2
C                            n1  n2
C
C  where N = n1+n2 and Z' means the transpose of Z. The first n1 columns
C  of Z span the specified invariant subspace of T.
C
C  If T has been obtained from the real Schur factorization of a matrix
C  A = Q*T*Q', then the reordered real Schur factorization of A is given
C  by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span
C  the corresponding invariant subspace of A.
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
      INTEGER           INFO, LDQ, LDT, LIWORK, LWORK, M, N
      CHARACTER         COMPQ, JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,*), T(LDT,*), WI(*), WORK(LWORK), WR(*)
      INTEGER           IWORK(LIWORK)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EST, RNORM, SCALE
      INTEGER           IERR, IFAIL, K, KASE, KK, KS, N1, N2, NN
      LOGICAL           PAIR, SWAP, WANTBH, WANTQ, WANTS, WANTSP
C     .. External Functions ..
      DOUBLE PRECISION  F06RAF
      EXTERNAL          F06RAF
C     .. External Subroutines ..
      EXTERNAL          DTREXC, DTRSYL, F04YCF, F06AAZ, F06QFF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
C     Decode and test the input parameters
C
      WANTBH = (JOB.EQ.'B' .OR. JOB.EQ.'b')
      WANTS = (JOB.EQ.'E' .OR. JOB.EQ.'e') .OR. WANTBH
      WANTSP = (JOB.EQ.'V' .OR. JOB.EQ.'v') .OR. WANTBH
      WANTQ = (COMPQ.EQ.'V' .OR. COMPQ.EQ.'v')
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
      ELSE
C
C        Set M to the dimension of the specified invariant subspace,
C        and test LWORK and LIWORK.
C
         M = 0
         PAIR = .FALSE.
         DO 20 K = 1, N
            IF (PAIR) THEN
               PAIR = .FALSE.
            ELSE
               IF (K.LT.N) THEN
                  IF (T(K+1,K).EQ.ZERO) THEN
                     IF (SELECT(K)) M = M + 1
                  ELSE
                     PAIR = .TRUE.
                     IF (SELECT(K) .OR. SELECT(K+1)) M = M + 2
                  END IF
               ELSE
                  IF (SELECT(N)) M = M + 1
               END IF
            END IF
   20    CONTINUE
C
         N1 = M
         N2 = N - M
         NN = N1*N2
C
         IF (LWORK.LT.1 .OR. ((WANTS .AND. .NOT. WANTSP)
     *       .AND. LWORK.LT.NN) .OR. (WANTSP .AND. LWORK.LT.2*NN)) THEN
            INFO = -15
         ELSE IF (LIWORK.LT.1 .OR. (WANTSP .AND. LIWORK.LT.NN)) THEN
            INFO = -17
         END IF
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QGF/DTRSEN',-INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF (M.EQ.N .OR. M.EQ.0) THEN
         IF (WANTS) S = ONE
         IF (WANTSP) SEP = F06RAF('1',N,N,T,LDT,WORK)
         GO TO 80
      END IF
C
C     Collect the selected blocks at the top-left corner of T.
C
      KS = 0
      PAIR = .FALSE.
      DO 40 K = 1, N
         IF (PAIR) THEN
            PAIR = .FALSE.
         ELSE
            SWAP = SELECT(K)
            IF (K.LT.N) THEN
               IF (T(K+1,K).NE.ZERO) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP .OR. SELECT(K+1)
               END IF
            END IF
            IF (SWAP) THEN
               KS = KS + 1
C
C              Swap the K-th block to position KS.
C
               IERR = 0
               KK = K
               IF (K.NE.KS) CALL DTREXC(COMPQ,N,T,LDT,Q,LDQ,KK,KS,WORK,
     *                                  IERR)
               IF (IERR.EQ.1 .OR. IERR.EQ.2) THEN
C
C                 Blocks too close to swap: exit.
C
                  INFO = 1
                  IF (WANTS) S = ZERO
                  IF (WANTSP) SEP = ZERO
                  GO TO 80
               END IF
               IF (PAIR) KS = KS + 1
            END IF
         END IF
   40 CONTINUE
C
      IF (WANTS) THEN
C
C        Solve Sylvester equation for R:
C
C           T11*R - R*T22 = scale*T12
C
         CALL F06QFF('General',N1,N2,T(1,N1+1),LDT,WORK,N1)
         CALL DTRSYL('N','N',-1,N1,N2,T,LDT,T(N1+1,N1+1),LDT,WORK,N1,
     *               SCALE,IERR)
C
C        Estimate the reciprocal of the condition number of the cluster
C        of eigenvalues.
C
         RNORM = F06RAF('F',N1,N2,WORK,N1,WORK)
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
         CALL F04YCF(KASE,NN,WORK,EST,WORK(NN+1),IWORK,IFAIL)
         IF (KASE.NE.0) THEN
            IF (KASE.EQ.1) THEN
C
C              Solve  T11*R - R*T22 = scale*X.
C
               CALL DTRSYL('N','N',-1,N1,N2,T,LDT,T(N1+1,N1+1),LDT,WORK,
     *                     N1,SCALE,IERR)
            ELSE
C
C              Solve  T11'*R - R*T22' = scale*X.
C
               CALL DTRSYL('T','T',-1,N1,N2,T,LDT,T(N1+1,N1+1),LDT,WORK,
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
C     Store the output eigenvalues in WR and WI.
C
      DO 100 K = 1, N
         WR(K) = T(K,K)
         WI(K) = ZERO
  100 CONTINUE
      DO 120 K = 1, N - 1
         IF (T(K+1,K).NE.ZERO) THEN
            WI(K) = SQRT(ABS(T(K,K+1)))*SQRT(ABS(T(K+1,K)))
            WI(K+1) = -WI(K)
         END IF
  120 CONTINUE
      RETURN
C
C     End of F08QGF (DTRSEN)
C
      END
