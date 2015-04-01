      SUBROUTINE F08QXF(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,MM,M,
     *                  WORK,RWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZTREVC(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,
     *                  LDVR,MM,M,WORK,RWORK,INFO)
C
C  Purpose
C  =======
C
C  ZTREVC computes all or some right and/or left eigenvectors of a
C  complex upper triangular matrix T.
C
C  The right eigenvector x and the left eigenvector y of T corresponding
C  to an eigenvalue w are defined by:
C
C               T*x = w*x,     y'*T = w*y'
C
C  where y' denotes the conjugate transpose of the vector y.
C
C  The routine may either return the matrices X and/or Y of right or
C  left eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
C  input unitary matrix. If T was obtained from the Schur factorization
C  of an original matrix A = Q*T*Q', then Q*X and/or Q*Y are the
C  matrices of right or left eigenvectors of A.
C
C  Arguments
C  =========
C
C  JOB     (input) CHARACTER*1
C          = 'R': compute right eigenvectors only
C          = 'L': compute left eigenvectors only
C          = 'B': compute both right and left eigenvectors
C
C  HOWMNY  (input) CHARACTER*1
C          Specifies how many left/right eigenvectors are wanted and the
C          form of the eigenvector matrix X or Y returned in VR or VL.
C          = 'A': compute all right and/or left eigenvectors;
C          = 'O': compute all right and/or left eigenvectors, multiplied
C                 on the left by an input (generally unitary) matrix;
C          = 'S': compute some right and/or left eigenvectors, specified
C                 by the logical array SELECT.
C
C  SELECT  (input) LOGICAL array, dimension (N)
C          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
C          computed.  To select the eigenvector corresponding to the
C          j-th eigenvalue, SELECT(j) must be set to .TRUE..
C          If HOWMNY = 'A' or 'O', SELECT is not referenced.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input/output) COMPLEX*16 array, dimension (LDT,N)
C          The upper triangular matrix T. T is modified, but restored.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)
C          On entry, if JOB = 'L' or 'B' and HOWMNY = 'O', VL must
C          contain an n-by-n matrix Q (usually the unitary matrix Q of
C          Schur vectors returned by ZHSEQR).
C          On exit, if JOB = 'L' or 'B', VL contains:
C          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
C          if HOWMNY = 'O', the matrix Q*Y;
C          if HOWMNY = 'S', the left eigenvectors of T specified by
C                           SELECT, stored consecutively in the columns
C                           of VL, in the same order as their
C                           eigenvalues.
C          If JOB = 'R', VL is not referenced.
C
C  LDVL    (input) INTEGER
C          The leading dimension of the array VL. LDVL >= max(1,N).
C
C  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
C          On entry, if JOB = 'R' or 'B' and HOWMNY = 'O', VR must
C          contain an n-by-n matrix Q (usually the unitary matrix Q of
C          Schur vectors returned by ZHSEQR).
C          On exit, if JOB = 'R' or 'B', VR contains:
C          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
C          if HOWMNY = 'O', the matrix Q*X;
C          if HOWMNY = 'S', the right eigenvectors of T specified by
C                           SELECT, stored consecutively in the columns
C                           of VR, in the same order as their
C                           eigenvalues.
C          If JOB = 'L', VR is not referenced.
C
C  LDVR    (input) INTEGER
C          The leading dimension of the array VR. LDVR >= max(1,N).
C
C  MM      (input) INTEGER
C          The number of columns in the arrays VL and/or VR. MM >= M.
C
C  M       (output) INTEGER
C          The number of columns in the arrays VL and/or VR required to
C          store the eigenvectors. If HOWMNY = 'A' or 'O', M is set
C          to N.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (2*N)
C
C  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:   successful exit
C          < 0:   if INFO = -i, the i-th argument had an illegal value
C
C  Further Details
C  ===============
C
C  The algorithm used in this program is basically backward (forward)
C  substitution, with scaling to make the the code robust against
C  possible overflow.
C
C  Each eigenvector is normalized so that the element of largest
C  magnitude has magnitude 1; here the magnitude of a complex number
C  (x,y) is taken to be |x| + |y|.
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
      COMPLEX*16        CMZERO, CMONE
      PARAMETER         (CMZERO=0.0D+0,CMONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDT, LDVL, LDVR, M, MM, N
      CHARACTER         HOWMNY, JOB
C     .. Array Arguments ..
      COMPLEX*16        T(LDT,*), VL(LDVL,*), VR(LDVR,*), WORK(*)
      DOUBLE PRECISION  RWORK(*)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM
      DOUBLE PRECISION  OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      INTEGER           I, II, IS, J, K, KI
      LOGICAL           ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV
C     .. External Functions ..
      DOUBLE PRECISION  DZASUM, X02AJF, X02AMF
      INTEGER           IZAMAX, X02BHF
      EXTERNAL          DZASUM, X02AJF, X02AMF, IZAMAX, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07TUZ, ZCOPY, ZDSCAL, ZGEMV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(CDUM) = ABS(DBLE(CDUM)) + ABS(DIMAG(CDUM))
C     .. Executable Statements ..
C
C     Decode and test the input parameters
C
      BOTHV = (JOB.EQ.'B' .OR. JOB.EQ.'b')
      RIGHTV = (JOB.EQ.'R' .OR. JOB.EQ.'r') .OR. BOTHV
      LEFTV = (JOB.EQ.'L' .OR. JOB.EQ.'l') .OR. BOTHV
C
      ALLV = (HOWMNY.EQ.'A' .OR. HOWMNY.EQ.'a')
      OVER = (HOWMNY.EQ.'O' .OR. HOWMNY.EQ.'o')
      SOMEV = (HOWMNY.EQ.'S' .OR. HOWMNY.EQ.'s')
C
C     Set M to the number of columns required to store the selected
C     eigenvectors.
C
      IF (SOMEV) THEN
         M = 0
         DO 20 J = 1, N
            IF (SELECT(J)) M = M + 1
   20    CONTINUE
      ELSE
         M = N
      END IF
C
      INFO = 0
      IF ( .NOT. RIGHTV .AND. .NOT. LEFTV) THEN
         INFO = -1
      ELSE IF ( .NOT. ALLV .AND. .NOT. OVER .AND. .NOT. SOMEV) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (LDT.LT.MAX(1,N)) THEN
         INFO = -6
      ELSE IF (LDVL.LT.1 .OR. (LEFTV .AND. LDVL.LT.N)) THEN
         INFO = -8
      ELSE IF (LDVR.LT.1 .OR. (RIGHTV .AND. LDVR.LT.N)) THEN
         INFO = -10
      ELSE IF (MM.LT.M) THEN
         INFO = -11
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QXF/ZTREVC',-INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF (N.EQ.0) RETURN
C
C     Set the constants to control overflow.
C
      UNFL = X02AMF()
      OVFL = ONE/UNFL
      ULP = X02AJF()*X02BHF()
      SMLNUM = UNFL*(N/ULP)
C
C     Store the diagonal elements of T in working array WORK.
C
      DO 40 I = 1, N
         WORK(I+N) = T(I,I)
   40 CONTINUE
C
C     Compute 1-norm of each column of strictly upper triangular
C     part of T to control overflow in triangular solver.
C
      RWORK(1) = ZERO
      DO 60 J = 2, N
         RWORK(J) = DZASUM(J-1,T(1,J),1)
   60 CONTINUE
C
      IF (RIGHTV) THEN
C
C        Compute right eigenvectors.
C
         IS = M
         DO 160 KI = N, 1, -1
C
            IF (SOMEV) THEN
               IF ( .NOT. SELECT(KI)) GO TO 160
            END IF
            SMIN = MAX(ULP*(CABS1(T(KI,KI))),SMLNUM)
C
            WORK(1) = CMONE
C
C           Form right-hand side.
C
            DO 80 K = 1, KI - 1
               WORK(K) = -T(K,KI)
   80       CONTINUE
C
C           Solve the triangular system:
C              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
C
            DO 100 K = 1, KI - 1
               T(K,K) = T(K,K) - T(KI,KI)
               IF (CABS1(T(K,K)).LT.SMIN) T(K,K) = SMIN
  100       CONTINUE
C
            IF (KI.GT.1) THEN
               CALL F07TUZ('Upper','No transpose','Non-unit','Y',KI-1,T,
     *                     LDT,WORK(1),SCALE,RWORK,INFO)
               WORK(KI) = SCALE
            END IF
C
C           Copy the vector x or Q*x to VR and normalize.
C
            IF ( .NOT. OVER) THEN
               CALL ZCOPY(KI,WORK(1),1,VR(1,IS),1)
C
               II = IZAMAX(KI,VR(1,IS),1)
               REMAX = ONE/CABS1(VR(II,IS))
               CALL ZDSCAL(KI,REMAX,VR(1,IS),1)
C
               DO 120 K = KI + 1, N
                  VR(K,IS) = CMZERO
  120          CONTINUE
            ELSE
               IF (KI.GT.1) CALL ZGEMV('N',N,KI-1,CMONE,VR,LDVR,WORK(1),
     *                                 1,DCMPLX(SCALE),VR(1,KI),1)
C
               II = IZAMAX(N,VR(1,KI),1)
               REMAX = ONE/CABS1(VR(II,KI))
               CALL ZDSCAL(N,REMAX,VR(1,KI),1)
            END IF
C
C           Set back the original diagonal elements of T.
C
            DO 140 K = 1, KI - 1
               T(K,K) = WORK(K+N)
  140       CONTINUE
C
            IS = IS - 1
  160    CONTINUE
      END IF
C
      IF (LEFTV) THEN
C
C        Compute left eigenvectors.
C
         IS = 1
         DO 260 KI = 1, N
C
            IF (SOMEV) THEN
               IF ( .NOT. SELECT(KI)) GO TO 260
            END IF
            SMIN = MAX(ULP*(CABS1(T(KI,KI))),SMLNUM)
C
            WORK(N) = CMONE
C
C           Form right-hand side.
C
            DO 180 K = KI + 1, N
               WORK(K) = -DCONJG(T(KI,K))
  180       CONTINUE
C
C           Solve the triangular system:
C              (T(KI+1:N,KI+1:N) - T(KI,KI))'*X = SCALE*WORK.
C
            DO 200 K = KI + 1, N
               T(K,K) = T(K,K) - T(KI,KI)
               IF (CABS1(T(K,K)).LT.SMIN) T(K,K) = SMIN
  200       CONTINUE
C
            IF (KI.LT.N) THEN
               CALL F07TUZ('Upper','Conjugate transpose','Non-unit','Y',
     *                     N-KI,T(KI+1,KI+1),LDT,WORK(KI+1),SCALE,RWORK,
     *                     INFO)
               WORK(KI) = SCALE
            END IF
C
C           Copy the vector x or Q*x to VL and normalize.
C
            IF ( .NOT. OVER) THEN
               CALL ZCOPY(N-KI+1,WORK(KI),1,VL(KI,IS),1)
C
               II = IZAMAX(N-KI+1,VL(KI,IS),1) + KI - 1
               REMAX = ONE/CABS1(VL(II,IS))
               CALL ZDSCAL(N-KI+1,REMAX,VL(KI,IS),1)
C
               DO 220 K = 1, KI - 1
                  VL(K,IS) = CMZERO
  220          CONTINUE
            ELSE
               IF (KI.LT.N) CALL ZGEMV('N',N,N-KI,CMONE,VL(1,KI+1),LDVL,
     *                                 WORK(KI+1),1,DCMPLX(SCALE),
     *                                 VL(1,KI),1)
C
               II = IZAMAX(N,VL(1,KI),1)
               REMAX = ONE/CABS1(VL(II,KI))
               CALL ZDSCAL(N,REMAX,VL(1,KI),1)
            END IF
C
C           Set back the original diagonal elements of T.
C
            DO 240 K = KI + 1, N
               T(K,K) = WORK(K+N)
  240       CONTINUE
C
            IS = IS + 1
  260    CONTINUE
      END IF
C
      RETURN
C
C     End of F08QXF (ZTREVC)
C
      END
