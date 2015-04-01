      SUBROUTINE F08PKF(JOB,EIGSRC,INITV,SELECT,N,H,LDH,WR,WI,VL,LDVL,
     *                  VR,LDVR,MM,M,WORK,IFAILL,IFAILR,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DHSEIN(JOB,EIGSRC,INITV,SELECT,N,H,LDH,WR,WI,VL,
     *                  LDVL,VR,LDVR,MM,M,WORK,IFAILL,IFAILR,INFO)
C
C  Purpose
C  =======
C
C  DHSEIN uses inverse iteration to find specified right and/or left
C  eigenvectors of a real upper Hessenberg matrix H.
C
C  The right eigenvector x and the left eigenvector y of the matrix H
C  corresponding to an eigenvalue w are defined by:
C
C               H x = w x,     y' H = w y'
C
C  where y' denotes the conjugate transpose of the vector y.
C
C  Arguments
C  =========
C
C  JOB     (input) CHARACTER*1
C          = 'R': compute right eigenvectors only;
C          = 'L': compute left eigenvectors only;
C          = 'B': compute both right and left eigenvectors.
C
C  EIGSRC  (input) CHARACTER*1
C          Specifies the source of eigenvalues supplied in WR and WI:
C          = 'Q': the eigenvalues were found using DHSEQR; thus, if
C                 H has zero sub-diagonal entries, and so is
C                 block-triangular, then the j-th eigenvalue can be
C                 assumed to be an eigenvalue of the block containing
C                 the j-th row/column.  This property allows DHSEIN to
C                 perform inverse iteration on just one diagonal block.
C          = 'N': no assumptions are made on the correspondence
C                 between eigenvalues and diagonal blocks.  In this
C                 case, DHSEIN must always perform inverse iteration
C                 using the whole matrix H.
C
C  INITV   (input) CHARACTER*1
C          Specifies whether initial starting vectors are supplied for
C          inverse iteration:
C          = 'N': no initial vectors are supplied;
C          = 'U': user-supplied initial vectors are stored in the arrays
C                 VL and/or VR.
C
C  SELECT  (input/output) LOGICAL array, dimension(N)
C          Specifies the eigenvectors to be computed. To select the
C          real eigenvector corresponding to a real eigenvalue WR(j),
C          SELECT(j) must be set to .TRUE.. To select the complex
C          eigenvector corresponding to a complex eigenvalue
C          (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)),
C          either SELECT(j) or SELECT(j+1) or both must be set to
C          .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is
C          .FALSE..
C
C  N       (input) INTEGER
C          The order of the matrix H.  N >= 0.
C
C  H       (input) DOUBLE PRECISION array, dimension (LDH,N)
C          The upper Hessenberg matrix H.
C
C  LDH     (input) INTEGER
C          The leading dimension of the array H.  LDH >= max(1,N).
C
C  WR      (input/output) DOUBLE PRECISION array, dimension (N)
C  WI      (input) DOUBLE PRECISION array, dimension (N)
C          On entry, the real and imaginary parts of the eigenvalues of
C          H; a complex conjugate pair of eigenvalues must be stored in
C          consecutive elements of WR and WI.
C          On exit, WR may have been altered since close eigenvalues
C          are perturbed slightly in searching for independent
C          eigenvectors.
C
C  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
C          On entry, if INITV = 'U' and JOB = 'L' or 'B', VL must
C          contain starting vectors for the inverse iteration for the
C          left eigenvectors; the starting vector for each eigenvector
C          must be in the same column(s) in which the eigenvector will
C          be stored.
C          On exit, if JOB = 'L' or 'B', the left eigenvectors
C          specified by SELECT will be stored consecutively in the
C          columns of VL, in the same order as their eigenvalues. A
C          complex eigenvector corresponding to a complex eigenvalue is
C          stored in two consecutive columns, the first holding the real
C          part and the second the imaginary part.
C          If JOB = 'R', VL is not referenced.
C
C  LDVL    (input) INTEGER
C          The leading dimension of the array VL.
C          LDVL >= max(1,N) if JOB = 'L' or 'B'; LDVL >= 1 otherwise.
C
C  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
C          On entry, if INITV = 'U' and JOB = 'R' or 'B', VR must
C          contain starting vectors for the inverse iteration for the
C          right eigenvectors; the starting vector for each eigenvector
C          must be in the same column(s) in which the eigenvector will
C          be stored.
C          On exit, if JOB = 'R' or 'B', the right eigenvectors
C          specified by SELECT will be stored consecutively in the
C          columns of VR, in the same order as their eigenvalues. A
C          complex eigenvector corresponding to a complex eigenvalue is
C          stored in two consecutive columns, the first holding the real
C          part and the second the imaginary part.
C          If JOB = 'L', VR is not referenced.
C
C  LDVR    (input) INTEGER
C          The leading dimension of the array VR.
C          LDVR >= max(1,N) if JOB = 'R' or 'B'; LDVR >= 1 otherwise.
C
C  MM      (input) INTEGER
C          The number of columns in the arrays VL and/or VR. MM >= M.
C
C  M       (output) INTEGER
C          The number of columns in the arrays VL and/or VR required to
C          store the eigenvectors; each selected real eigenvector
C          occupies one column and each selected complex eigenvector
C          occupies two columns.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension ((N+2)*N)
C
C  IFAILL  (output) INTEGER array, dimension (MM)
C          If JOB = 'L' or 'B', IFAILL(i) = j > 0 if the left
C          eigenvector in the i-th column of VL (corresponding to the
C          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the
C          eigenvector converged satisfactorily. If the i-th and (i+1)th
C          columns of VL hold a complex eigenvector, then IFAILL(i) and
C          IFAILL(i+1) are set to the same value.
C          If JOB = 'R', IFAILL is not referenced.
C
C  IFAILR  (output) INTEGER array, dimension (MM)
C          If JOB = 'R' or 'B', IFAILR(i) = j > 0 if the right
C          eigenvector in the i-th column of VR (corresponding to the
C          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the
C          eigenvector converged satisfactorily. If the i-th and (i+1)th
C          columns of VR hold a complex eigenvector, then IFAILR(i) and
C          IFAILR(i+1) are set to the same value.
C          If JOB = 'L', IFAILR is not referenced.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C          > 0: INFO is the number of eigenvectors which failed to
C               converge; see IFAILL and IFAILR for further details.
C
C  Further Details
C  ===============
C
C  Each eigenvector is normalized so that the element of largest
C  magnitude has magnitude 1; here the magnitude of a complex number
C  (x,y) is taken to be |x|+|y|.
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
      INTEGER           INFO, LDH, LDVL, LDVR, M, MM, N
      CHARACTER         EIGSRC, INITV, JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LDH,*), VL(LDVL,*), VR(LDVR,*), WI(*),
     *                  WORK(*), WR(*)
      INTEGER           IFAILL(*), IFAILR(*)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGNUM, EPS3, HNORM, SMLNUM, ULP, UNFL, WKI, WKR
      INTEGER           I, IINFO, K, KL, KLN, KR, KSI, KSR, LDWORK
      LOGICAL           BOTHV, FROMQR, LEFTV, NOINIT, PAIR, RIGHTV
C     .. External Functions ..
      DOUBLE PRECISION  F06RMF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          F06RMF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08PKZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
C     Decode and test the input parameters.
C
      BOTHV = (JOB.EQ.'B' .OR. JOB.EQ.'b')
      RIGHTV = (JOB.EQ.'R' .OR. JOB.EQ.'r') .OR. BOTHV
      LEFTV = (JOB.EQ.'L' .OR. JOB.EQ.'l') .OR. BOTHV
C
      FROMQR = (EIGSRC.EQ.'Q' .OR. EIGSRC.EQ.'q')
C
      NOINIT = (INITV.EQ.'N' .OR. INITV.EQ.'n')
C
C     Set M to the number of columns required to store the selected
C     eigenvectors, and standardize the array SELECT.
C
      M = 0
      PAIR = .FALSE.
      DO 20 K = 1, N
         IF (PAIR) THEN
            PAIR = .FALSE.
            SELECT(K) = .FALSE.
         ELSE
            IF (WI(K).EQ.ZERO) THEN
               IF (SELECT(K)) M = M + 1
            ELSE
               PAIR = .TRUE.
               IF (SELECT(K) .OR. SELECT(K+1)) THEN
                  SELECT(K) = .TRUE.
                  M = M + 2
               END IF
            END IF
         END IF
   20 CONTINUE
C
      INFO = 0
      IF ( .NOT. RIGHTV .AND. .NOT. LEFTV) THEN
         INFO = -1
      ELSE IF ( .NOT. FROMQR .AND. .NOT.
     *         (EIGSRC.EQ.'N' .OR. EIGSRC.EQ.'n')) THEN
         INFO = -2
      ELSE IF ( .NOT. NOINIT .AND. .NOT.
     *         (INITV.EQ.'U' .OR. INITV.EQ.'u')) THEN
         INFO = -3
      ELSE IF (N.LT.0) THEN
         INFO = -5
      ELSE IF (LDH.LT.MAX(1,N)) THEN
         INFO = -7
      ELSE IF (LDVL.LT.1 .OR. (LEFTV .AND. LDVL.LT.N)) THEN
         INFO = -11
      ELSE IF (LDVR.LT.1 .OR. (RIGHTV .AND. LDVR.LT.N)) THEN
         INFO = -13
      ELSE IF (MM.LT.M) THEN
         INFO = -14
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08PKF/DHSEIN',-INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF (N.EQ.0) RETURN
C
C     Set machine-dependent constants.
C
      UNFL = X02AMF()
      ULP = X02AJF()*X02BHF()
      SMLNUM = UNFL*(N/ULP)
      BIGNUM = (ONE-ULP)/SMLNUM
C
      LDWORK = N + 1
C
      KL = 1
      KLN = 0
      IF (FROMQR) THEN
         KR = 0
      ELSE
         KR = N
      END IF
      KSR = 1
C
      DO 240 K = 1, N
         IF (SELECT(K)) THEN
C
C           Compute eigenvector(s) corresponding to W(K).
C
            IF (FROMQR) THEN
C
C              If affiliation of eigenvalues is known, check whether
C              the matrix splits.
C
C              Determine KL and KR such that 1 <= KL <= K <= KR <= N
C              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
C              KR = N).
C
C              Then inverse iteration can be performed with the
C              submatrix H(KL:N,KL:N) for a left eigenvector, and with
C              the submatrix H(1:KR,1:KR) for a right eigenvector.
C
               DO 40 I = K, KL + 1, -1
                  IF (H(I,I-1).EQ.ZERO) GO TO 60
   40          CONTINUE
   60          CONTINUE
               KL = I
               IF (K.GT.KR) THEN
                  DO 80 I = K, N - 1
                     IF (H(I+1,I).EQ.ZERO) GO TO 100
   80             CONTINUE
  100             CONTINUE
                  KR = I
               END IF
            END IF
C
            IF (KL.NE.KLN) THEN
               KLN = KL
C
C              Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
C              has not ben computed before.
C
               HNORM = F06RMF('I',KR-KL+1,H(KL,KL),LDH,WORK)
               IF (HNORM.GT.ZERO) THEN
                  EPS3 = HNORM*ULP
               ELSE
                  EPS3 = SMLNUM
               END IF
            END IF
C
C           Perturb eigenvalue if it is close to any previous
C           selected eigenvalues affiliated to the submatrix
C           H(KL:KR,KL:KR). Close roots are modified by EPS3.
C
            WKR = WR(K)
            WKI = WI(K)
  120       CONTINUE
            DO 140 I = K - 1, KL, -1
               IF (SELECT(I) .AND. ABS(WR(I)-WKR)+ABS(WI(I)-WKI)
     *             .LT.EPS3) THEN
                  WKR = WKR + EPS3
                  GO TO 120
               END IF
  140       CONTINUE
            WR(K) = WKR
C
            PAIR = WKI .NE. ZERO
            IF (PAIR) THEN
               KSI = KSR + 1
            ELSE
               KSI = KSR
            END IF
            IF (LEFTV) THEN
C
C              Compute left eigenvector.
C
               CALL F08PKZ(.FALSE.,NOINIT,N-KL+1,H(KL,KL),LDH,WKR,WKI,
     *                     VL(KL,KSR),VL(KL,KSI),WORK,LDWORK,
     *                     WORK(N*N+N+1),EPS3,SMLNUM,BIGNUM,IINFO)
               IF (IINFO.GT.0) THEN
                  IF (PAIR) THEN
                     INFO = INFO + 2
                  ELSE
                     INFO = INFO + 1
                  END IF
                  IFAILL(KSR) = K
                  IFAILL(KSI) = K
               ELSE
                  IFAILL(KSR) = 0
                  IFAILL(KSI) = 0
               END IF
               DO 160 I = 1, KL - 1
                  VL(I,KSR) = ZERO
  160          CONTINUE
               IF (PAIR) THEN
                  DO 180 I = 1, KL - 1
                     VL(I,KSI) = ZERO
  180             CONTINUE
               END IF
            END IF
            IF (RIGHTV) THEN
C
C              Compute right eigenvector.
C
               CALL F08PKZ(.TRUE.,NOINIT,KR,H,LDH,WKR,WKI,VR(1,KSR),
     *                     VR(1,KSI),WORK,LDWORK,WORK(N*N+N+1),EPS3,
     *                     SMLNUM,BIGNUM,IINFO)
               IF (IINFO.GT.0) THEN
                  IF (PAIR) THEN
                     INFO = INFO + 2
                  ELSE
                     INFO = INFO + 1
                  END IF
                  IFAILR(KSR) = K
                  IFAILR(KSI) = K
               ELSE
                  IFAILR(KSR) = 0
                  IFAILR(KSI) = 0
               END IF
               DO 200 I = KR + 1, N
                  VR(I,KSR) = ZERO
  200          CONTINUE
               IF (PAIR) THEN
                  DO 220 I = KR + 1, N
                     VR(I,KSI) = ZERO
  220             CONTINUE
               END IF
            END IF
C
            IF (PAIR) THEN
               KSR = KSR + 2
            ELSE
               KSR = KSR + 1
            END IF
         END IF
  240 CONTINUE
C
      RETURN
C
C     End of F08PKF (DHSEIN)
C
      END
