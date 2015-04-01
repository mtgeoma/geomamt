      SUBROUTINE F08QYF(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,S,SEP,
     *                  MM,M,WORK,LDWORK,RWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZTRSNA(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,
     *                  LDVR,S,SEP,MM,M,WORK,LDWORK,RWORK,INFO)
C
C  Purpose
C  =======
C
C  ZTRSNA estimates reciprocal condition numbers for specified
C  eigenvalues and/or right eigenvectors of a complex upper triangular
C  matrix T (or of any matrix Q*T*Q' with Q unitary).
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
C          for the j-th eigenpair, SELECT(j) must be set to .TRUE..
C          If HOWMNY = 'A', SELECT is not referenced.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input) COMPLEX*16 array, dimension (LDT,N)
C          The upper triangular matrix T.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  VL      (input) COMPLEX*16 array, dimension (LDVL,M)
C          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
C          (or of any Q*T*Q' with Q unitary), corresponding to the
C          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
C          must be stored in consecutive columns of VL, as returned by
C          ZHSEIN or ZTREVC.
C          If JOB = 'V', VL is not referenced.
C
C  LDVL    (input) INTEGER
C          The leading dimension of the array VL.
C          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
C
C  VR      (input) COMPLEX*16 array, dimension (LDVR,M)
C          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
C          (or of any Q*T*Q' with Q unitary), corresponding to the
C          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
C          must be stored in consecutive columns of VR, as returned by
C          ZHSEIN or ZTREVC.
C          If JOB = 'V', VR is not referenced.
C
C  LDVR    (input) INTEGER
C          The leading dimension of the array VR.
C          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
C
C  S       (output) DOUBLE PRECISION array, dimension (MM)
C          If JOB = 'E' or 'B', the reciprocal condition numbers of the
C          selected eigenvalues, stored in consecutive elements of the
C          array. Thus S(j), SEP(j), and the j-th columns of VL and VR
C          all correspond to the same eigenpair (but not in general the
C          j-th eigenpair, unless all eigenpairs are selected).
C          If JOB = 'V', S is not referenced.
C
C  SEP     (output) DOUBLE PRECISION array, dimension (MM)
C          If JOB = 'V' or 'B', the estimated reciprocal condition
C          numbers of the selected eigenvectors, stored in consecutive
C          elements of the array.
C          If JOB = 'E', SEP is not referenced.
C
C  MM      (input) INTEGER
C          The number of elements in the arrays S and SEP. MM >= M.
C
C  M       (output) INTEGER
C          The number of elements of the arrays S and SEP used to store
C          the specified condition numbers. If HOWMNY = 'A', M is set
C          to N.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,N+1)
C          If JOB = 'E', WORK is not referenced.
C
C  LDWORK  (input) INTEGER
C          The leading dimension of the array WORK.
C          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
C
C  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
C          If JOB = 'E', RWORK is not referenced.
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
C  to lambda; v' denotes the conjugate transpose of v, and norm(u)
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
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D0+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
      CHARACTER         HOWMNY, JOB
C     .. Array Arguments ..
      COMPLEX*16        T(LDT,*), VL(LDVL,*), VR(LDVR,*), WORK(LDWORK,*)
      DOUBLE PRECISION  RWORK(*), S(*), SEP(*)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM, PROD
      DOUBLE PRECISION  BIGNUM, EPS, EST, LNRM, RNRM, SCALE, SMLNUM,
     *                  XNORM
      INTEGER           I, IERR, IFAIL, IX, J, K, KASE, KS
      LOGICAL           SOMCON, WANTBH, WANTS, WANTSP
      CHARACTER         NORMIN
C     .. Local Arrays ..
      COMPLEX*16        DUMMY(1)
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      DOUBLE PRECISION  DZNRM2, X02AJF, X02AMF
      INTEGER           IZAMAX, X02BHF
      EXTERNAL          ZDOTC, DZNRM2, X02AJF, X02AMF, IZAMAX, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F04ZCF, F06AAZ, F06TFF, F07AUZ, F07TUZ, ZTREXC
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DIMAG, MAX
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(CDUM) = ABS(DBLE(CDUM)) + ABS(DIMAG(CDUM))
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
C     Set M to the number of eigenpairs for which condition numbers are
C     to be computed.
C
      IF (SOMCON) THEN
         M = 0
         DO 20 J = 1, N
            IF (SELECT(J)) M = M + 1
   20    CONTINUE
      ELSE
         M = N
      END IF
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
      ELSE IF (MM.LT.M) THEN
         INFO = -13
      ELSE IF (LDWORK.LT.1 .OR. (WANTSP .AND. LDWORK.LT.N)) THEN
         INFO = -16
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QYF/ZTRSNA',-INFO)
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
      KS = 1
      DO 100 K = 1, N
C
         IF (SOMCON) THEN
            IF ( .NOT. SELECT(K)) GO TO 100
         END IF
C
         IF (WANTS) THEN
C
C           Compute the reciprocal condition number of the k-th
C           eigenvalue.
C
            PROD = ZDOTC(N,VR(1,KS),1,VL(1,KS),1)
            RNRM = DZNRM2(N,VR(1,KS),1)
            LNRM = DZNRM2(N,VL(1,KS),1)
            S(KS) = ABS(PROD)/(RNRM*LNRM)
C
         END IF
C
         IF (WANTSP) THEN
C
C           Estimate the reciprocal condition number of the k-th
C           eigenvector.
C
C           Copy the matrix T to the array WORK and swap the k-th
C           diagonal element to the (1,1) position.
C
            CALL F06TFF('General',N,N,T,LDT,WORK,LDWORK)
            CALL ZTREXC('No Q',N,WORK,LDWORK,DUMMY,1,K,1,IERR)
C
C           Form  C = T22 - lambda*I in WORK(2:N,2:N).
C
            DO 40 I = 2, N
               WORK(I,I) = WORK(I,I) - WORK(1,1)
   40       CONTINUE
C
C           Estimate a lower bound for the 1-norm of inv(C'). The 1st
C           and (N+1)th columns of WORK are used to store work vectors.
C
            SEP(KS) = ZERO
            EST = ZERO
            KASE = 0
            NORMIN = 'N'
   60       CONTINUE
            IFAIL = 0
            CALL F04ZCF(KASE,N-1,WORK,EST,WORK(1,N+1),IFAIL)
C
            IF (KASE.NE.0) THEN
               IF (KASE.EQ.1) THEN
C
C                 Solve C'*x = scale*b
C
                  CALL F07TUZ('Upper','Conjugate transpose','Nonunit',
     *                        NORMIN,N-1,WORK(2,2),LDWORK,WORK,SCALE,
     *                        RWORK,IERR)
               ELSE
C
C                 Solve C*x = scale*b
C
                  CALL F07TUZ('Upper','No transpose','Nonunit',NORMIN,
     *                        N-1,WORK(2,2),LDWORK,WORK,SCALE,RWORK,
     *                        IERR)
               END IF
               NORMIN = 'Y'
               IF (SCALE.NE.ONE) THEN
C
C                 Multiply by 1/SCALE if doing so will not cause
C                 overflow.
C
                  IX = IZAMAX(N-1,WORK,1)
                  XNORM = CABS1(WORK(IX,1))
                  IF (SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO)
     *                GO TO 80
                  CALL F07AUZ(N,SCALE,WORK,1)
               END IF
               GO TO 60
            END IF
C
            SEP(KS) = ONE/MAX(EST,SMLNUM)
         END IF
C
   80    CONTINUE
         KS = KS + 1
  100 CONTINUE
      RETURN
C
C     End of F08QYF (ZTRSNA)
C
      END
