      SUBROUTINE F08PXF(JOB,EIGSRC,INITV,SELECT,N,H,LDH,W,VL,LDVL,VR,
     *                  LDVR,MM,M,WORK,RWORK,IFAILL,IFAILR,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZHSEIN(JOB,EIGSRC,INITV,SELECT,N,H,LDH,W,VL,
     *                  LDVL,VR,LDVR,MM,M,WORK,RWORK,IFAILL,IFAILR,INFO)
C
C  Purpose
C  =======
C
C  ZHSEIN uses inverse iteration to find specified right and/or left
C  eigenvectors of a complex upper Hessenberg matrix H.
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
C          Specifies the source of eigenvalues supplied in W:
C          = 'Q': the eigenvalues were found using ZHSEQR; thus, if
C                 H has zero sub-diagonal entries, and so is
C                 block-triangular, then the j-th eigenvalue can be
C                 assumed to be an eigenvalue of the block containing
C                 the j-th row/column.  This property allows ZHSEIN to
C                 perform inverse iteration on just one diagonal block.
C          = 'N': no assumptions are made on the correspondence
C                 between eigenvalues and diagonal blocks.  In this
C                 case, ZHSEIN must always perform inverse iteration
C                 using the whole matrix H.
C
C  INITV   (input) CHARACTER*1
C          Specifies whether initial starting vectors are supplied for
C          inverse iteration:
C          = 'N': no initial vectors are supplied;
C          = 'U': user-supplied initial vectors are stored in the arrays
C                 VL and/or VR.
C
C  SELECT  (input) LOGICAL array, dimension (N)
C          Specifies the eigenvectors to be computed. To select the
C          eigenvector corresponding to the eigenvalue W(j),
C          SELECT(j) must be set to .TRUE..
C
C  N       (input) INTEGER
C          The order of the matrix H.  N >= 0.
C
C  H       (input) COMPLEX*16 array, dimension (LDH,N)
C          The upper Hessenberg matrix H.
C
C  LDH     (input) INTEGER
C          The leading dimension of the array H.  LDH >= max(1,N).
C
C  W       (input/output) COMPLEX*16 array, dimension (N)
C          On entry, the eigenvalues of H.
C          On exit, the real parts of W may have been altered since
C          close eigenvalues are perturbed slightly in searching for
C          independent eigenvectors.
C
C  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)
C          On entry, if INITV = 'U' and JOB = 'L' or 'B', VL must
C          contain starting vectors for the inverse iteration for the
C          left eigenvectors; the starting vector for each eigenvector
C          must be in the same column in which the eigenvector will be
C          stored.
C          On exit, if JOB = 'L' or 'B', the left eigenvectors
C          specified by SELECT will be stored consecutively in the
C          columns of VL, in the same order as their eigenvalues.
C          If JOB = 'R', VL is not referenced.
C
C  LDVL    (input) INTEGER
C          The leading dimension of the array VL.
C          LDVL >= max(1,N) if JOB = 'L' or 'B'; LDVL >= 1 otherwise.
C
C  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
C          On entry, if INITV = 'U' and JOB = 'R' or 'B', VR must
C          contain starting vectors for the inverse iteration for the
C          right eigenvectors; the starting vector for each eigenvector
C          must be in the same column in which the eigenvector will be
C          stored.
C          On exit, if JOB = 'R' or 'B', the right eigenvectors
C          specified by SELECT will be stored consecutively in the
C          columns of VR, in the same order as their eigenvalues.
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
C          store the eigenvectors (= the number of .TRUE. elements in
C          SELECT).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N*N)
C
C  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
C
C  IFAILL  (output) INTEGER array, dimension (MM)
C          If JOB = 'L' or 'B', IFAILL(i) = j > 0 if the left
C          eigenvector in the i-th column of VL (corresponding to the
C          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the
C          eigenvector converged satisfactorily.
C          If JOB = 'R', IFAILL is not referenced.
C
C  IFAILR  (output) INTEGER array, dimension (MM)
C          If JOB = 'R' or 'B', IFAILR(i) = j > 0 if the right
C          eigenvector in the i-th column of VR (corresponding to the
C          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the
C          eigenvector converged satisfactorily.
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
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  RZERO
      PARAMETER         (RZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDH, LDVL, LDVR, M, MM, N
      CHARACTER         EIGSRC, INITV, JOB
C     .. Array Arguments ..
      COMPLEX*16        H(LDH,*), VL(LDVL,*), VR(LDVR,*), W(*), WORK(*)
      DOUBLE PRECISION  RWORK(*)
      INTEGER           IFAILL(*), IFAILR(*)
      LOGICAL           SELECT(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM, WK
      DOUBLE PRECISION  EPS3, HNORM, SMLNUM, ULP, UNFL
      INTEGER           I, IINFO, K, KL, KLN, KR, KS, LDWORK
      LOGICAL           BOTHV, FROMQR, LEFTV, NOINIT, RIGHTV
C     .. External Functions ..
      DOUBLE PRECISION  F06UMF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          F06UMF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08PXZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DIMAG, MAX
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(CDUM) = ABS(DBLE(CDUM)) + ABS(DIMAG(CDUM))
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
C     eigenvectors.
C
      M = 0
      DO 20 K = 1, N
         IF (SELECT(K)) M = M + 1
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
         INFO = -10
      ELSE IF (LDVR.LT.1 .OR. (RIGHTV .AND. LDVR.LT.N)) THEN
         INFO = -12
      ELSE IF (MM.LT.M) THEN
         INFO = -13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08PXF/ZHSEIN',-INFO)
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
C
      LDWORK = N
C
      KL = 1
      KLN = 0
      IF (FROMQR) THEN
         KR = 0
      ELSE
         KR = N
      END IF
      KS = 1
C
      DO 200 K = 1, N
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
               HNORM = F06UMF('I',KR-KL+1,H(KL,KL),LDH,RWORK)
               IF (HNORM.GT.RZERO) THEN
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
            WK = W(K)
  120       CONTINUE
            DO 140 I = K - 1, KL, -1
               IF (SELECT(I) .AND. CABS1(W(I)-WK).LT.EPS3) THEN
                  WK = WK + EPS3
                  GO TO 120
               END IF
  140       CONTINUE
            W(K) = WK
C
            IF (LEFTV) THEN
C
C              Compute left eigenvector.
C
               CALL F08PXZ(.FALSE.,NOINIT,N-KL+1,H(KL,KL),LDH,WK,
     *                     VL(KL,KS),WORK,LDWORK,RWORK,EPS3,SMLNUM,
     *                     IINFO)
               IF (IINFO.GT.0) THEN
                  INFO = INFO + 1
                  IFAILL(KS) = K
               ELSE
                  IFAILL(KS) = 0
               END IF
               DO 160 I = 1, KL - 1
                  VL(I,KS) = ZERO
  160          CONTINUE
            END IF
            IF (RIGHTV) THEN
C
C              Compute right eigenvector.
C
               CALL F08PXZ(.TRUE.,NOINIT,KR,H,LDH,WK,VR(1,KS),WORK,
     *                     LDWORK,RWORK,EPS3,SMLNUM,IINFO)
               IF (IINFO.GT.0) THEN
                  INFO = INFO + 1
                  IFAILR(KS) = K
               ELSE
                  IFAILR(KS) = 0
               END IF
               DO 180 I = KR + 1, N
                  VR(I,KS) = ZERO
  180          CONTINUE
            END IF
            KS = KS + 1
         END IF
  200 CONTINUE
C
      RETURN
C
C     End of F08PXF (ZHSEIN)
C
      END
