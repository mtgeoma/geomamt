      SUBROUTINE F08QKF(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,MM,M,
     *                  WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DTREVC(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,
     *                  LDVR,MM,M,WORK,INFO)
C
C  Purpose
C  =======
C
C  DTREVC computes all or some right and/or left eigenvectors of a
C  real upper quasi-triangular matrix T in Schur canonical form.
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
C  input orthogonal matrix. If T was obtained from the real Schur
C  factorization of an original matrix A = Q*T*Q', then Q*X and/or Q*Y
C  are the matrices of right or left eigenvectors of A.
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
C                 on the left by an input (generally orthogonal) matrix;
C          = 'S': compute some right and/or left eigenvectors, specified
C                 by the logical array SELECT.
C
C  SELECT  (input/output) LOGICAL array, dimension (N)
C          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
C          computed. To select the real eigenvector corresponding to a
C          real eigenvalue w(j), SELECT(j) must be set to .TRUE.. To
C          select the complex eigenvector corresponding to a complex
C          conjugate pair w(j) and w(j+1), either SELECT(j) or
C          SELECT(j+1) must be set to .TRUE.; then on exit SELECT(j) is
C          .TRUE. and SELECT(j+1) is .FALSE..
C          If HOWMNY = 'A' or 'O', SELECT is not referenced.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
C          The upper quasi-triangular matrix T in Schur canonical form.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
C          On entry, if JOB = 'L' or 'B' and HOWMNY = 'O', VL must
C          contain an n-by-n matrix Q (usually the orthogonal matrix Q
C          of Schur vectors returned by DHSEQR).
C          On exit, if JOB = 'L' or 'B', VL contains:
C          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
C          if HOWMNY = 'O', the matrix Q*Y;
C          if HOWMNY = 'S', the left eigenvectors of T specified by
C                           SELECT, stored consecutively in the columns
C                           of VL, in the same order as their
C                           eigenvalues.
C          A complex eigenvector corresponding to a complex eigenvalue
C          is stored in two consecutive columns, the first holding the
C          real part, and the second the imaginary part.
C          If JOB = 'R', VL is not referenced.
C
C  LDVL    (input) INTEGER
C          The leading dimension of the array VL. LDVL >= max(1,N).
C
C  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
C          On entry, if JOB = 'R' or 'B' and HOWMNY = 'O', VR must
C          contain an n-by-n matrix Q (usually the orthogonal matrix Q
C          of Schur vectors returned by DHSEQR).
C          On exit, if JOB = 'R' or 'B', VR contains:
C          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
C          if HOWMNY = 'O', the matrix Q*X;
C          if HOWMNY = 'S', the right eigenvectors of T specified by
C                           SELECT, stored consecutively in the columns
C                           of VR, in the same order as their
C                           eigenvalues.
C          A complex eigenvector corresponding to a complex eigenvalue
C          is stored in two consecutive columns, the first holding the
C          real part and the second the imaginary part.
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
C          store the eigenvectors; each selected real eigenvector
C          occupies one column and each selected complex eigenvector
C          occupies two columns.  If HOWMNY = 'A' or 'O', M is set to N.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
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
C     .. Scalar Arguments ..
      INTEGER           INFO, LDT, LDVL, LDVR, M, MM, N
      CHARACTER         HOWMNY, JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  T(LDT,*), VL(LDVL,*), VR(LDVR,*), WORK(*)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE,
     *                  SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR,
     *                  XNORM
      INTEGER           I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, N2
      LOGICAL           ALLV, BOTHV, LEFTV, OVER, PAIR, RIGHTV, SOMEV
C     .. Local Arrays ..
      DOUBLE PRECISION  X(2,2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF, X02AMF
      INTEGER           IDAMAX, X02BHF
      EXTERNAL          DDOT, X02AJF, X02AMF, IDAMAX, X02BHF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, F06AAZ, F08QHX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
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
      ELSE
C
C        Set M to the number of columns required to store the selected
C        eigenvectors, standardize the array SELECT if necessary, and
C        test MM.
C
         IF (SOMEV) THEN
            M = 0
            PAIR = .FALSE.
            DO 20 J = 1, N
               IF (PAIR) THEN
                  PAIR = .FALSE.
                  SELECT(J) = .FALSE.
               ELSE
                  IF (J.LT.N) THEN
                     IF (T(J+1,J).EQ.ZERO) THEN
                        IF (SELECT(J)) M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF (SELECT(J) .OR. SELECT(J+1)) THEN
                           SELECT(J) = .TRUE.
                           M = M + 2
                        END IF
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
            INFO = -11
         END IF
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QKF/DTREVC',-INFO)
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
      BIGNUM = (ONE-ULP)/SMLNUM
C
C     Compute 1-norm of each column of strictly upper triangular
C     part of T to control overflow in triangular solver.
C
      WORK(1) = ZERO
      DO 60 J = 2, N
         WORK(J) = ZERO
         DO 40 I = 1, J - 1
            WORK(J) = WORK(J) + ABS(T(I,J))
   40    CONTINUE
   60 CONTINUE
C
C     Index IP is used to specify the real or complex eigenvalue:
C       IP = 0, real eigenvalue,
C            1, first of conjugate complex pair: (wr,wi)
C           -1, second of conjugate complex pair: (wr,wi)
C
      N2 = 2*N
C
      IF (RIGHTV) THEN
C
C        Compute right eigenvectors.
C
         IP = 0
         IS = M
         DO 280 KI = N, 1, -1
C
            IF (IP.EQ.1) GO TO 260
            IF (KI.EQ.1) GO TO 80
            IF (T(KI,KI-1).EQ.ZERO) GO TO 80
            IP = -1
C
   80       CONTINUE
            IF (SOMEV) THEN
               IF (IP.EQ.0) THEN
                  IF ( .NOT. SELECT(KI)) GO TO 260
               ELSE
                  IF ( .NOT. SELECT(KI-1)) GO TO 260
               END IF
            END IF
C
C           Compute the KI-th eigenvalue (WR,WI).
C
            WR = T(KI,KI)
            WI = ZERO
            IF (IP.NE.0) WI = SQRT(ABS(T(KI,KI-1)))*SQRT(ABS(T(KI-1,KI))
     *                        )
            SMIN = MAX(ULP*(ABS(WR)+ABS(WI)),SMLNUM)
C
            IF (IP.EQ.0) THEN
C
C              Real right eigenvector
C
               WORK(KI+N) = ONE
C
C              Form right-hand side
C
               DO 100 K = 1, KI - 1
                  WORK(K+N) = -T(K,KI)
  100          CONTINUE
C
C              Solve the upper quasi-triangular system:
C                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
C
               JNXT = KI - 1
               DO 120 J = KI - 1, 1, -1
                  IF (J.GT.JNXT) GO TO 120
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF (J.GT.1) THEN
                     IF (T(J,J-1).NE.ZERO) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
C
                  IF (J1.EQ.J2) THEN
C
C                    1-by-1 diagonal block
C
                     CALL F08QHX(.FALSE.,1,1,SMIN,ONE,T(J,J),LDT,ONE,
     *                           ONE,WORK(J+N),N,WR,ZERO,X,2,SCALE,
     *                           XNORM,IERR)
C
C                    Scale X(1,1) to avoid overflow when updating
C                    the right-hand side.
C
                     IF (XNORM.GT.ONE) THEN
                        IF (WORK(J).GT.BIGNUM/XNORM) THEN
                           X(1,1) = X(1,1)/XNORM
                           SCALE = SCALE/XNORM
                        END IF
                     END IF
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) CALL DSCAL(KI,SCALE,WORK(1+N),1)
                     WORK(J+N) = X(1,1)
C
C                    Update right-hand side
C
                     CALL DAXPY(J-1,-X(1,1),T(1,J),1,WORK(1+N),1)
C
                  ELSE
C
C                    2-by-2 diagonal block
C
                     CALL F08QHX(.FALSE.,2,1,SMIN,ONE,T(J-1,J-1),LDT,
     *                           ONE,ONE,WORK(J-1+N),N,WR,ZERO,X,2,
     *                           SCALE,XNORM,IERR)
C
C                    Scale X(1,1) and X(2,1) to avoid overflow when
C                    updating the right-hand side.
C
                     IF (XNORM.GT.ONE) THEN
                        BETA = MAX(WORK(J-1),WORK(J))
                        IF (BETA.GT.BIGNUM/XNORM) THEN
                           X(1,1) = X(1,1)/XNORM
                           X(2,1) = X(2,1)/XNORM
                           SCALE = SCALE/XNORM
                        END IF
                     END IF
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) CALL DSCAL(KI,SCALE,WORK(1+N),1)
                     WORK(J-1+N) = X(1,1)
                     WORK(J+N) = X(2,1)
C
C                    Update right-hand side
C
                     CALL DAXPY(J-2,-X(1,1),T(1,J-1),1,WORK(1+N),1)
                     CALL DAXPY(J-2,-X(2,1),T(1,J),1,WORK(1+N),1)
                  END IF
  120          CONTINUE
C
C              Copy the vector x or Q*x to VR and normalize.
C
               IF ( .NOT. OVER) THEN
                  CALL DCOPY(KI,WORK(1+N),1,VR(1,IS),1)
C
                  II = IDAMAX(KI,VR(1,IS),1)
                  REMAX = ONE/ABS(VR(II,IS))
                  CALL DSCAL(KI,REMAX,VR(1,IS),1)
C
                  DO 140 K = KI + 1, N
                     VR(K,IS) = ZERO
  140             CONTINUE
               ELSE
                  IF (KI.GT.1) CALL DGEMV('N',N,KI-1,ONE,VR,LDVR,
     *                                    WORK(1+N),1,WORK(KI+N),
     *                                    VR(1,KI),1)
C
                  II = IDAMAX(N,VR(1,KI),1)
                  REMAX = ONE/ABS(VR(II,KI))
                  CALL DSCAL(N,REMAX,VR(1,KI),1)
               END IF
C
            ELSE
C
C              Complex right eigenvector.
C
C              Initial solve
C                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
C                [ (T(KI,KI-1)   T(KI,KI)   )               ]
C
               IF (ABS(T(KI-1,KI)).GE.ABS(T(KI,KI-1))) THEN
                  WORK(KI-1+N) = ONE
                  WORK(KI+N2) = WI/T(KI-1,KI)
               ELSE
                  WORK(KI-1+N) = -WI/T(KI,KI-1)
                  WORK(KI+N2) = ONE
               END IF
               WORK(KI+N) = ZERO
               WORK(KI-1+N2) = ZERO
C
C              Form right-hand side
C
               DO 160 K = 1, KI - 2
                  WORK(K+N) = -WORK(KI-1+N)*T(K,KI-1)
                  WORK(K+N2) = -WORK(KI+N2)*T(K,KI)
  160          CONTINUE
C
C              Solve upper quasi-triangular system:
C              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
C
               JNXT = KI - 2
               DO 180 J = KI - 2, 1, -1
                  IF (J.GT.JNXT) GO TO 180
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF (J.GT.1) THEN
                     IF (T(J,J-1).NE.ZERO) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
C
                  IF (J1.EQ.J2) THEN
C
C                    1-by-1 diagonal block
C
                     CALL F08QHX(.FALSE.,1,2,SMIN,ONE,T(J,J),LDT,ONE,
     *                           ONE,WORK(J+N),N,WR,WI,X,2,SCALE,XNORM,
     *                           IERR)
C
C                    Scale X(1,1) and X(1,2) to avoid overflow when
C                    updating the right-hand side.
C
                     IF (XNORM.GT.ONE) THEN
                        IF (WORK(J).GT.BIGNUM/XNORM) THEN
                           X(1,1) = X(1,1)/XNORM
                           X(1,2) = X(1,2)/XNORM
                           SCALE = SCALE/XNORM
                        END IF
                     END IF
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) THEN
                        CALL DSCAL(KI,SCALE,WORK(1+N),1)
                        CALL DSCAL(KI,SCALE,WORK(1+N2),1)
                     END IF
                     WORK(J+N) = X(1,1)
                     WORK(J+N2) = X(1,2)
C
C                    Update the right-hand side
C
                     CALL DAXPY(J-1,-X(1,1),T(1,J),1,WORK(1+N),1)
                     CALL DAXPY(J-1,-X(1,2),T(1,J),1,WORK(1+N2),1)
C
                  ELSE
C
C                    2-by-2 diagonal block
C
                     CALL F08QHX(.FALSE.,2,2,SMIN,ONE,T(J-1,J-1),LDT,
     *                           ONE,ONE,WORK(J-1+N),N,WR,WI,X,2,SCALE,
     *                           XNORM,IERR)
C
C                    Scale X to avoid overflow when updating
C                    the right-hand side.
C
                     IF (XNORM.GT.ONE) THEN
                        BETA = MAX(WORK(J-1),WORK(J))
                        IF (BETA.GT.BIGNUM/XNORM) THEN
                           REC = ONE/XNORM
                           X(1,1) = X(1,1)*REC
                           X(1,2) = X(1,2)*REC
                           X(2,1) = X(2,1)*REC
                           X(2,2) = X(2,2)*REC
                           SCALE = SCALE*REC
                        END IF
                     END IF
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) THEN
                        CALL DSCAL(KI,SCALE,WORK(1+N),1)
                        CALL DSCAL(KI,SCALE,WORK(1+N2),1)
                     END IF
                     WORK(J-1+N) = X(1,1)
                     WORK(J+N) = X(2,1)
                     WORK(J-1+N2) = X(1,2)
                     WORK(J+N2) = X(2,2)
C
C                    Update the right-hand side
C
                     CALL DAXPY(J-2,-X(1,1),T(1,J-1),1,WORK(1+N),1)
                     CALL DAXPY(J-2,-X(2,1),T(1,J),1,WORK(1+N),1)
                     CALL DAXPY(J-2,-X(1,2),T(1,J-1),1,WORK(1+N2),1)
                     CALL DAXPY(J-2,-X(2,2),T(1,J),1,WORK(1+N2),1)
                  END IF
  180          CONTINUE
C
C              Copy the vector x or Q*x to VR and normalize.
C
               IF ( .NOT. OVER) THEN
                  CALL DCOPY(KI,WORK(1+N),1,VR(1,IS-1),1)
                  CALL DCOPY(KI,WORK(1+N2),1,VR(1,IS),1)
C
                  EMAX = ZERO
                  DO 200 K = 1, KI
                     EMAX = MAX(EMAX,ABS(VR(K,IS-1))+ABS(VR(K,IS)))
  200             CONTINUE
C
                  REMAX = ONE/EMAX
                  CALL DSCAL(KI,REMAX,VR(1,IS-1),1)
                  CALL DSCAL(KI,REMAX,VR(1,IS),1)
C
                  DO 220 K = KI + 1, N
                     VR(K,IS-1) = ZERO
                     VR(K,IS) = ZERO
  220             CONTINUE
C
               ELSE
C
                  IF (KI.GT.2) THEN
                     CALL DGEMV('N',N,KI-2,ONE,VR,LDVR,WORK(1+N),1,
     *                          WORK(KI-1+N),VR(1,KI-1),1)
                     CALL DGEMV('N',N,KI-2,ONE,VR,LDVR,WORK(1+N2),1,
     *                          WORK(KI+N2),VR(1,KI),1)
                  ELSE
                     CALL DSCAL(N,WORK(KI-1+N),VR(1,KI-1),1)
                     CALL DSCAL(N,WORK(KI+N2),VR(1,KI),1)
                  END IF
C
                  EMAX = ZERO
                  DO 240 K = 1, N
                     EMAX = MAX(EMAX,ABS(VR(K,KI-1))+ABS(VR(K,KI)))
  240             CONTINUE
                  REMAX = ONE/EMAX
                  CALL DSCAL(N,REMAX,VR(1,KI-1),1)
                  CALL DSCAL(N,REMAX,VR(1,KI),1)
               END IF
            END IF
C
            IS = IS - 1
            IF (IP.NE.0) IS = IS - 1
  260       CONTINUE
            IF (IP.EQ.1) IP = 0
            IF (IP.EQ.-1) IP = 1
  280    CONTINUE
      END IF
C
      IF (LEFTV) THEN
C
C        Compute left eigenvectors.
C
         IP = 0
         IS = 1
         DO 500 KI = 1, N
C
            IF (IP.EQ.-1) GO TO 480
            IF (KI.EQ.N) GO TO 300
            IF (T(KI+1,KI).EQ.ZERO) GO TO 300
            IP = 1
C
  300       CONTINUE
            IF (SOMEV) THEN
               IF ( .NOT. SELECT(KI)) GO TO 480
            END IF
C
C           Compute the KI-th eigenvalue (WR,WI).
C
            WR = T(KI,KI)
            WI = ZERO
            IF (IP.NE.0) WI = SQRT(ABS(T(KI,KI+1)))*SQRT(ABS(T(KI+1,KI))
     *                        )
            SMIN = MAX(ULP*(ABS(WR)+ABS(WI)),SMLNUM)
C
            IF (IP.EQ.0) THEN
C
C              Real left eigenvector.
C
               WORK(KI+N) = ONE
C
C              Form right-hand side
C
               DO 320 K = KI + 1, N
                  WORK(K+N) = -T(KI,K)
  320          CONTINUE
C
C              Solve the quasi-triangular system:
C                 (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
C
               VMAX = ONE
               VCRIT = BIGNUM
C
               JNXT = KI + 1
               DO 340 J = KI + 1, N
                  IF (J.LT.JNXT) GO TO 340
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF (J.LT.N) THEN
                     IF (T(J+1,J).NE.ZERO) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
C
                  IF (J1.EQ.J2) THEN
C
C                    1-by-1 diagonal block
C
C                    Scale if necessary to avoid overflow when forming
C                    the right-hand side.
C
                     IF (WORK(J).GT.VCRIT) THEN
                        REC = ONE/VMAX
                        CALL DSCAL(N-KI+1,REC,WORK(KI+N),1)
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
C
                     WORK(J+N) = WORK(J+N) - DDOT(J-KI-1,T(KI+1,J),1,
     *                           WORK(KI+1+N),1)
C
C                    Solve (T(J,J)-WR)'*X = WORK
C
                     CALL F08QHX(.FALSE.,1,1,SMIN,ONE,T(J,J),LDT,ONE,
     *                           ONE,WORK(J+N),N,WR,ZERO,X,2,SCALE,
     *                           XNORM,IERR)
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) CALL DSCAL(N-KI+1,SCALE,
     *                                      WORK(KI+N),1)
                     WORK(J+N) = X(1,1)
                     VMAX = MAX(ABS(WORK(J+N)),VMAX)
                     VCRIT = BIGNUM/VMAX
C
                  ELSE
C
C                    2-by-2 diagonal block
C
C                    Scale if necessary to avoid overflow when forming
C                    the right-hand side.
C
                     BETA = MAX(WORK(J),WORK(J+1))
                     IF (BETA.GT.VCRIT) THEN
                        REC = ONE/VMAX
                        CALL DSCAL(N-KI+1,REC,WORK(KI+N),1)
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
C
                     WORK(J+N) = WORK(J+N) - DDOT(J-KI-1,T(KI+1,J),1,
     *                           WORK(KI+1+N),1)
C
                     WORK(J+1+N) = WORK(J+1+N) - DDOT(J-KI-1,T(KI+1,J+1)
     *                             ,1,WORK(KI+1+N),1)
C
C                    Solve
C                      [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
C                      [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
C
                     CALL F08QHX(.TRUE.,2,1,SMIN,ONE,T(J,J),LDT,ONE,ONE,
     *                           WORK(J+N),N,WR,ZERO,X,2,SCALE,XNORM,
     *                           IERR)
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) CALL DSCAL(N-KI+1,SCALE,
     *                                      WORK(KI+N),1)
                     WORK(J+N) = X(1,1)
                     WORK(J+1+N) = X(2,1)
C
                     VMAX = MAX(ABS(WORK(J+N)),ABS(WORK(J+1+N)),VMAX)
                     VCRIT = BIGNUM/VMAX
C
                  END IF
  340          CONTINUE
C
C              Copy the vector x or Q*x to VL and normalize.
C
               IF ( .NOT. OVER) THEN
                  CALL DCOPY(N-KI+1,WORK(KI+N),1,VL(KI,IS),1)
C
                  II = IDAMAX(N-KI+1,VL(KI,IS),1) + KI - 1
                  REMAX = ONE/ABS(VL(II,IS))
                  CALL DSCAL(N-KI+1,REMAX,VL(KI,IS),1)
C
                  DO 360 K = 1, KI - 1
                     VL(K,IS) = ZERO
  360             CONTINUE
C
               ELSE
C
                  IF (KI.LT.N) CALL DGEMV('N',N,N-KI,ONE,VL(1,KI+1),
     *                                    LDVL,WORK(KI+1+N),1,WORK(KI+N)
     *                                    ,VL(1,KI),1)
C
                  II = IDAMAX(N,VL(1,KI),1)
                  REMAX = ONE/ABS(VL(II,KI))
                  CALL DSCAL(N,REMAX,VL(1,KI),1)
C
               END IF
C
            ELSE
C
C              Complex left eigenvector.
C
C               Initial solve:
C                 ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
C                 ((T(KI+1,KI) T(KI+1,KI+1))                )
C
               IF (ABS(T(KI,KI+1)).GE.ABS(T(KI+1,KI))) THEN
                  WORK(KI+N) = WI/T(KI,KI+1)
                  WORK(KI+1+N2) = ONE
               ELSE
                  WORK(KI+N) = ONE
                  WORK(KI+1+N2) = -WI/T(KI+1,KI)
               END IF
               WORK(KI+1+N) = ZERO
               WORK(KI+N2) = ZERO
C
C              Form right-hand side
C
               DO 380 K = KI + 2, N
                  WORK(K+N) = -WORK(KI+N)*T(KI,K)
                  WORK(K+N2) = -WORK(KI+1+N2)*T(KI+1,K)
  380          CONTINUE
C
C              Solve complex quasi-triangular system:
C              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
C
               VMAX = ONE
               VCRIT = BIGNUM
C
               JNXT = KI + 2
               DO 400 J = KI + 2, N
                  IF (J.LT.JNXT) GO TO 400
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF (J.LT.N) THEN
                     IF (T(J+1,J).NE.ZERO) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
C
                  IF (J1.EQ.J2) THEN
C
C                    1-by-1 diagonal block
C
C                    Scale if necessary to avoid overflow when
C                    forming the right-hand side entries.
C
                     IF (WORK(J).GT.VCRIT) THEN
                        REC = ONE/VMAX
                        CALL DSCAL(N-KI+1,REC,WORK(KI+N),1)
                        CALL DSCAL(N-KI+1,REC,WORK(KI+N2),1)
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
C
                     WORK(J+N) = WORK(J+N) - DDOT(J-KI-2,T(KI+2,J),1,
     *                           WORK(KI+2+N),1)
                     WORK(J+N2) = WORK(J+N2) - DDOT(J-KI-2,T(KI+2,J),1,
     *                            WORK(KI+2+N2),1)
C
C                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
C
                     CALL F08QHX(.FALSE.,1,2,SMIN,ONE,T(J,J),LDT,ONE,
     *                           ONE,WORK(J+N),N,WR,-WI,X,2,SCALE,XNORM,
     *                           IERR)
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) THEN
                        CALL DSCAL(N-KI+1,SCALE,WORK(KI+N),1)
                        CALL DSCAL(N-KI+1,SCALE,WORK(KI+N2),1)
                     END IF
                     WORK(J+N) = X(1,1)
                     WORK(J+N2) = X(1,2)
                     VMAX = MAX(ABS(WORK(J+N)),ABS(WORK(J+N2)),VMAX)
                     VCRIT = BIGNUM/VMAX
C
                  ELSE
C
C                    2-by-2 diagonal block
C
C                    Scale if necessary to avoid overflow when forming
C                    the right-hand side entries.
C
                     BETA = MAX(WORK(J),WORK(J+1))
                     IF (BETA.GT.VCRIT) THEN
                        REC = ONE/VMAX
                        CALL DSCAL(N-KI+1,REC,WORK(KI+N),1)
                        CALL DSCAL(N-KI+1,REC,WORK(KI+N2),1)
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
C
                     WORK(J+N) = WORK(J+N) - DDOT(J-KI-2,T(KI+2,J),1,
     *                           WORK(KI+2+N),1)
C
                     WORK(J+N2) = WORK(J+N2) - DDOT(J-KI-2,T(KI+2,J),1,
     *                            WORK(KI+2+N2),1)
C
                     WORK(J+1+N) = WORK(J+1+N) - DDOT(J-KI-2,T(KI+2,J+1)
     *                             ,1,WORK(KI+2+N),1)
C
                     WORK(J+1+N2) = WORK(J+1+N2) - DDOT(J-KI-2,
     *                              T(KI+2,J+1),1,WORK(KI+2+N2),1)
C
C                    Solve 2-by-2 complex linear equation
C                      ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
C                      ([T(j+1,j) T(j+1,j+1)]             )
C
                     CALL F08QHX(.TRUE.,2,2,SMIN,ONE,T(J,J),LDT,ONE,ONE,
     *                           WORK(J+N),N,WR,-WI,X,2,SCALE,XNORM,
     *                           IERR)
C
C                    Scale if necessary
C
                     IF (SCALE.NE.ONE) THEN
                        CALL DSCAL(N-KI+1,SCALE,WORK(KI+N),1)
                        CALL DSCAL(N-KI+1,SCALE,WORK(KI+N2),1)
                     END IF
                     WORK(J+N) = X(1,1)
                     WORK(J+N2) = X(1,2)
                     WORK(J+1+N) = X(2,1)
                     WORK(J+1+N2) = X(2,2)
                     VMAX = MAX(ABS(X(1,1)),ABS(X(1,2)),ABS(X(2,1)),
     *                      ABS(X(2,2)),VMAX)
                     VCRIT = BIGNUM/VMAX
C
                  END IF
  400          CONTINUE
C
C              Copy the vector x or Q*x to VL and normalize.
C
               IF ( .NOT. OVER) THEN
                  CALL DCOPY(N-KI+1,WORK(KI+N),1,VL(KI,IS),1)
                  CALL DCOPY(N-KI+1,WORK(KI+N2),1,VL(KI,IS+1),1)
C
                  EMAX = ZERO
                  DO 420 K = KI, N
                     EMAX = MAX(EMAX,ABS(VL(K,IS))+ABS(VL(K,IS+1)))
  420             CONTINUE
                  REMAX = ONE/EMAX
                  CALL DSCAL(N-KI+1,REMAX,VL(KI,IS),1)
                  CALL DSCAL(N-KI+1,REMAX,VL(KI,IS+1),1)
C
                  DO 440 K = 1, KI - 1
                     VL(K,IS) = ZERO
                     VL(K,IS+1) = ZERO
  440             CONTINUE
               ELSE
                  IF (KI.LT.N-1) THEN
                     CALL DGEMV('N',N,N-KI-1,ONE,VL(1,KI+2),LDVL,
     *                          WORK(KI+2+N),1,WORK(KI+N),VL(1,KI),1)
                     CALL DGEMV('N',N,N-KI-1,ONE,VL(1,KI+2),LDVL,
     *                          WORK(KI+2+N2),1,WORK(KI+1+N2),VL(1,KI+1)
     *                          ,1)
                  ELSE
                     CALL DSCAL(N,WORK(KI+N),VL(1,KI),1)
                     CALL DSCAL(N,WORK(KI+1+N2),VL(1,KI+1),1)
                  END IF
C
                  EMAX = ZERO
                  DO 460 K = 1, N
                     EMAX = MAX(EMAX,ABS(VL(K,KI))+ABS(VL(K,KI+1)))
  460             CONTINUE
                  REMAX = ONE/EMAX
                  CALL DSCAL(N,REMAX,VL(1,KI),1)
                  CALL DSCAL(N,REMAX,VL(1,KI+1),1)
C
               END IF
C
            END IF
C
            IS = IS + 1
            IF (IP.NE.0) IS = IS + 1
  480       CONTINUE
            IF (IP.EQ.-1) IP = 0
            IF (IP.EQ.1) IP = -1
C
  500    CONTINUE
C
      END IF
C
      RETURN
C
C     End of F08QKF (DTREVC)
C
      END
