      SUBROUTINE F08QLF(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,S,SEP,
     *                  MM,M,WORK,LDWORK,IWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DTRSNA(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,
     *                  LDVR,S,SEP,MM,M,WORK,LDWORK,IWORK,INFO)
C
C  Purpose
C  =======
C
C  DTRSNA estimates reciprocal condition numbers for specified
C  eigenvalues and/or right eigenvectors of a real upper
C  quasi-triangular matrix T (or of any matrix Q*T*Q' with Q
C  orthogonal).
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
C          Specifies whether condition numbers are required for
C          eigenvalues (S) or eigenvectors (SEP):
C          = 'E': for eigenvalues only (S);
C          = 'V': for eigenvectors only (SEP);
C          = 'B': for both eigenvalues and eigenvectors (S and SEP).
C
C  HOWMNY  (input) CHARACTER*1
C          = 'A': compute condition numbers for all eigenpairs;
C          = 'S': compute condition numbers for selected eigenpairs
C                 specified by the array SELECT.
C
C  SELECT  (input) LOGICAL array, dimension (N)
C          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
C          condition numbers are required. To select condition numbers
C          for the eigenpair corresponding to a real eigenvalue w(j),
C          SELECT(j) must be set to .TRUE.. To select condition numbers
C          corresponding to a complex conjugate pair of eigenvalues w(j)
C          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
C          set to .TRUE..
C          If HOWMNY = 'A', SELECT is not referenced.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
C          The upper quasi-triangular matrix T, in Schur canonical form.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  VL      (input) DOUBLE PRECISION array, dimension (LDVL,M)
C          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
C          (or of any Q*T*Q' with Q orthogonal), corresponding to the
C          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
C          must be stored in consecutive columns of VL, as returned by
C          DHSEIN or DTREVC.
C          If JOB = 'V', VL is not referenced.
C
C  LDVL    (input) INTEGER
C          The leading dimension of the array VL.
C          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
C
C  VR      (input) DOUBLE PRECISION array, dimension (LDVR,M)
C          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
C          (or of any Q*T*Q' with Q orthogonal), corresponding to the
C          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
C          must be stored in consecutive columns of VR, as returned by
C          DHSEIN or DTREVC.
C          If JOB = 'V', VR is not referenced.
C
C  LDVR    (input) INTEGER
C          The leading dimension of the array VR.
C          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
C
C  S       (output) DOUBLE PRECISION array, dimension (MM)
C          If JOB = 'E' or 'B', the reciprocal condition numbers of the
C          selected eigenvalues, stored in consecutive elements of the
C          array. For a complex conjugate pair of eigenvalues two
C          consecutive elements of S are set to the same value. Thus
C          S(j), SEP(j), and the j-th columns of VL and VR all
C          correspond to the same eigenpair (but not in general the
C          j-th eigenpair, unless all eigenpairs are selected).
C          If JOB = 'V', S is not referenced.
C
C  SEP     (output) DOUBLE PRECISION array, dimension (MM)
C          If JOB = 'V' or 'B', the estimated reciprocal condition
C          numbers of the selected eigenvectors, stored in consecutive
C          elements of the array. For a complex eigenvector two
C          consecutive elements of SEP are set to the same value. If
C          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)
C          is set to 0; this can only occur when the true value would be
C          very small anyway.
C          If JOB = 'E', SEP is not referenced.
C
C  MM      (input) INTEGER
C          The number of elements in the arrays S and SEP. MM >= M.
C
C  M       (output) INTEGER
C          The number of elements of the arrays S and SEP used to store
C          the specified condition numbers; for each selected real
C          eigenvalue one element is used, and for each selected complex
C          conjugate pair of eigenvalues, two elements are used. If
C          HOWMNY = 'A', M is set to N.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,N+1)
C          If JOB = 'E', WORK is not referenced.
C
C  LDWORK  (input) INTEGER
C          The leading dimension of the array WORK.
C          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
C
C  IWORK   (workspace) INTEGER array, dimension (N)
C          If JOB = 'E', IWORK is not referenced.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  Further Details
C  ===============
C
C  The reciprocal of the condition number of an eigenvalue lambda is
C  defined as
C
C          S(lambda) = |v'*u| / (norm(u)*norm(v))
C
C  where u and v are the right and left eigenvectors of T corresponding
C  to lambda; v' denotes the conjugate-transpose of v, and norm(u)
C  denotes the Euclidean norm. These reciprocal condition numbers always
C  lie between zero (very badly conditioned) and one (very well
C  conditioned). If n = 1, S(lambda) is defined to be 1.
C
C  An approximate error bound for a computed eigenvalue W(i) is given by
C
C                      EPS * norm(T) / S(i)
C
C  where EPS is the machine precision.
C
C  The reciprocal of the condition number of the right eigenvector u
C  corresponding to lambda is defined as follows. Suppose
C
C              T = ( lambda  c  )
C                  (   0    T22 )
C
C  Then the reciprocal condition number is
C
C          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
C
C  where sigma-min denotes the smallest singular value. We approximate
C  the smallest singular value by the reciprocal of an estimate of the
C  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
C  defined to be abs(T(1,1)).
C
C  An approximate error bound for a computed right eigenvector VR(i)
C  is given by
C
C                      EPS * norm(T) / SEP(i)
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
      CHARACTER         HOWMNY, JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  S(*), SEP(*), T(LDT,*), VL(LDVL,*), VR(LDVR,*),
     *                  WORK(LDWORK,*)
      INTEGER           IWORK(*)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGNUM, COND, CS, DELTA, DUMM, EPS, EST, LNRM,
     *                  MU, PROD, PROD1, PROD2, RNRM, SCALE, SMLNUM, SN
      INTEGER           I, IERR, IFAIL, IFST, ILST, J, K, KASE, KS, N2,
     *                  NN
      LOGICAL           PAIR, SOMCON, WANTBH, WANTS, WANTSP
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BNF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          DDOT, DNRM2, F06BNF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          DTREXC, F04YCF, F06AAZ, F06QFF, F08QLZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
C     Decode and test the input parameters
C
      WANTBH = (JOB.EQ.'B' .OR. JOB.EQ.'b')
      WANTS = (JOB.EQ.'E' .OR. JOB.EQ.'e') .OR. WANTBH
      WANTSP = (JOB.EQ.'V' .OR. JOB.EQ.'v') .OR. WANTBH
C
      SOMCON = (HOWMNY.EQ.'S' .OR. HOWMNY.EQ.'s')
C
      INFO = 0
      IF ( .NOT. WANTS .AND. .NOT. WANTSP) THEN
         INFO = -1
      ELSE IF ( .NOT. (HOWMNY.EQ.'A' .OR. HOWMNY.EQ.'a')
     *         .AND. .NOT. SOMCON) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (LDT.LT.MAX(1,N)) THEN
         INFO = -6
      ELSE IF (LDVL.LT.1 .OR. (WANTS .AND. LDVL.LT.N)) THEN
         INFO = -8
      ELSE IF (LDVR.LT.1 .OR. (WANTS .AND. LDVR.LT.N)) THEN
         INFO = -10
      ELSE
C
C        Set M to the number of eigenpairs for which condition numbers
C        are required, and test MM.
C
         IF (SOMCON) THEN
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
   20       CONTINUE
         ELSE
            M = N
         END IF
C
         IF (MM.LT.M) THEN
            INFO = -13
         ELSE IF (LDWORK.LT.1 .OR. (WANTSP .AND. LDWORK.LT.N)) THEN
            INFO = -16
         END IF
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QLF/DTRSNA',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF (N.EQ.1) THEN
         IF (SOMCON) THEN
            IF ( .NOT. SELECT(1)) RETURN
         END IF
         IF (WANTS) S(1) = ONE
         IF (WANTSP) SEP(1) = ABS(T(1,1))
         RETURN
      END IF
C
C     Get machine constants
C
      EPS = X02AJF()*X02BHF()
      SMLNUM = X02AMF()/EPS
      BIGNUM = ONE/SMLNUM
C
      KS = 0
      PAIR = .FALSE.
      DO 120 K = 1, N
C
C        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.
C
         IF (PAIR) THEN
            PAIR = .FALSE.
            GO TO 120
         ELSE
            IF (K.LT.N) PAIR = T(K+1,K) .NE. ZERO
         END IF
C
C        Determine whether condition numbers are required for the k-th
C        eigenpair.
C
         IF (SOMCON) THEN
            IF (PAIR) THEN
               IF ( .NOT. SELECT(K) .AND. .NOT. SELECT(K+1)) GO TO 120
            ELSE
               IF ( .NOT. SELECT(K)) GO TO 120
            END IF
         END IF
C
         KS = KS + 1
C
         IF (WANTS) THEN
C
C           Compute the reciprocal condition number of the k-th
C           eigenvalue.
C
            IF ( .NOT. PAIR) THEN
C
C              Real eigenvalue.
C
               PROD = DDOT(N,VR(1,KS),1,VL(1,KS),1)
               RNRM = DNRM2(N,VR(1,KS),1)
               LNRM = DNRM2(N,VL(1,KS),1)
               S(KS) = ABS(PROD)/(RNRM*LNRM)
            ELSE
C
C              Complex eigenvalue.
C
               PROD1 = DDOT(N,VR(1,KS),1,VL(1,KS),1)
               PROD1 = PROD1 + DDOT(N,VR(1,KS+1),1,VL(1,KS+1),1)
               PROD2 = DDOT(N,VL(1,KS),1,VR(1,KS+1),1)
               PROD2 = PROD2 - DDOT(N,VL(1,KS+1),1,VR(1,KS),1)
               RNRM = F06BNF(DNRM2(N,VR(1,KS),1),DNRM2(N,VR(1,KS+1),1))
               LNRM = F06BNF(DNRM2(N,VL(1,KS),1),DNRM2(N,VL(1,KS+1),1))
               COND = F06BNF(PROD1,PROD2)/(RNRM*LNRM)
               S(KS) = COND
               S(KS+1) = COND
            END IF
         END IF
C
         IF (WANTSP) THEN
C
C           Estimate the reciprocal condition number of the k-th
C           eigenvector.
C
C           Copy the matrix T to the array WORK and swap the diagonal
C           block beginning at T(k,k) to the (1,1) position.
C
            CALL F06QFF('General',N,N,T,LDT,WORK,LDWORK)
            IFST = K
            ILST = 1
            CALL DTREXC('No Q',N,WORK,LDWORK,DUMMY,1,IFST,ILST,
     *                  WORK(1,N+1),IERR)
C
            IF (IERR.EQ.1 .OR. IERR.EQ.2) THEN
C
C              Could not swap because blocks not well separated
C
               SCALE = ONE
               EST = BIGNUM
            ELSE
C
C              Reordering successful
C
               IF (WORK(2,1).EQ.ZERO) THEN
C
C                 Form C = T22 - lambda*I in WORK(2:N,2:N).
C
                  DO 40 I = 2, N
                     WORK(I,I) = WORK(I,I) - WORK(1,1)
   40             CONTINUE
                  N2 = 1
                  NN = N - 1
               ELSE
C
C                 Triangularize the 2 by 2 block by unitary
C                 transformation U = [  cs   i*ss ]
C                                    [ i*ss   cs  ].
C                 such that the (1,1) position of WORK is complex
C                 eigenvalue lambda with positive imaginary part. (2,2)
C                 position of WORK is the complex eigenvalue lambda
C                 with negative imaginary  part.
C
                  MU = SQRT(ABS(WORK(1,2)))*SQRT(ABS(WORK(2,1)))
                  DELTA = F06BNF(MU,WORK(2,1))
                  CS = MU/DELTA
                  SN = -WORK(2,1)/DELTA
C
C                 Form
C
C                 C' = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
C                                        [   mu                     ]
C                                        [         ..               ]
C                                        [             ..           ]
C                                        [                  mu      ]
C                 where C' is conjugate transpose of complex matrix C,
C                 and RWORK is stored starting in the N+1-st column of
C                 WORK.
C
                  DO 60 J = 3, N
                     WORK(2,J) = CS*WORK(2,J)
                     WORK(J,J) = WORK(J,J) - WORK(1,1)
   60             CONTINUE
                  WORK(2,2) = ZERO
C
                  WORK(1,N+1) = TWO*MU
                  DO 80 I = 2, N - 1
                     WORK(I,N+1) = SN*WORK(1,I+1)
   80             CONTINUE
                  N2 = 2
                  NN = 2*(N-1)
               END IF
C
C              Estimate norm(inv(C'))
C
               EST = ZERO
               KASE = 0
  100          CONTINUE
               IFAIL = 0
               CALL F04YCF(KASE,NN,WORK(1,N+4),EST,WORK(1,N+2),IWORK,
     *                     IFAIL)
               IF (KASE.NE.0) THEN
                  IF (KASE.EQ.1) THEN
                     IF (N2.EQ.1) THEN
C
C                       Real eigenvalue: solve C'*x = scale*c.
C
                        CALL F08QLZ(.TRUE.,.TRUE.,N-1,WORK(2,2),LDWORK,
     *                              DUMMY,DUMM,SCALE,WORK(1,N+4),
     *                              WORK(1,N+6),IERR)
                     ELSE
C
C                       Complex eigenvalue: solve
C                       C'*(p+iq) = scale*(c+id) in real arithmetic.
C
                        CALL F08QLZ(.TRUE.,.FALSE.,N-1,WORK(2,2),LDWORK,
     *                              WORK(1,N+1),MU,SCALE,WORK(1,N+4),
     *                              WORK(1,N+6),IERR)
                     END IF
                  ELSE
                     IF (N2.EQ.1) THEN
C
C                       Real eigenvalue: solve C*x = scale*c.
C
                        CALL F08QLZ(.FALSE.,.TRUE.,N-1,WORK(2,2),LDWORK,
     *                              DUMMY,DUMM,SCALE,WORK(1,N+4),
     *                              WORK(1,N+6),IERR)
                     ELSE
C
C                       Complex eigenvalue: solve
C                       C*(p+iq) = scale*(c+id) in real arithmetic.
C
                        CALL F08QLZ(.FALSE.,.FALSE.,N-1,WORK(2,2),
     *                              LDWORK,WORK(1,N+1),MU,SCALE,
     *                              WORK(1,N+4),WORK(1,N+6),IERR)
C
                     END IF
                  END IF
C
                  GO TO 100
               END IF
            END IF
C
            SEP(KS) = SCALE/MAX(EST,SMLNUM)
            IF (PAIR) SEP(KS+1) = SEP(KS)
         END IF
C
         IF (PAIR) KS = KS + 1
C
  120 CONTINUE
      RETURN
C
C     End of F08QLF (DTRSNA)
C
      END
