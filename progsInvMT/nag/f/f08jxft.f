      SUBROUTINE F08JXF(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAIL,
     *                  INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1655 (JUN 1995).
C     .. Entry Points ..
      ENTRY             ZSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,
     *                  IFAIL,INFO)
C
C  Purpose
C  =======
C
C  ZSTEIN computes the eigenvectors of a real symmetric tridiagonal
C  matrix corresponding to specified eigenvalues, using inverse
C  iteration.  
C
C  The maximum number of iterations allowed for each eigenvector is
C  specified by an internal parameter MAXITS (currently set to 5).
C
C  Although the eigenvectors are real, they are stored in a complex
C  array, which may be passed to ZUNMTR or ZUPMTR for back
C  transformation to the eigenvectors of a complex Hermitian matrix
C  which was reduced to tridiagonal form.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix.  N >= 0.
C
C  D       (input) DOUBLE PRECISION array, dimension (N)
C          The n diagonal elements of the tridiagonal matrix T.
C
C  E       (input) DOUBLE PRECISION array, dimension (N)
C          The (n-1) subdiagonal elements of the tridiagonal matrix
C          T, in elements 1 to N-1.  E(N) need not be set.
C
C  M       (input) INTEGER
C          The number of eigenvectors to be found.  0 <= M <= N.
C
C  W       (input) DOUBLE PRECISION array, dimension (N)
C          The first M elements of W contain the eigenvalues for
C          which eigenvectors are to be computed.  The eigenvalues
C          should be grouped by split-off block and ordered from
C          smallest to largest within the block.  ( The output array
C          W from DSTEBZ with ORDER = 'B' is expected here. )
C
C  IBLOCK  (input) INTEGER array, dimension (N)
C          The submatrix indices associated with the corresponding
C          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
C          the first submatrix from the top, =2 if W(i) belongs to
C          the second submatrix, etc.  ( The output array IBLOCK
C          from DSTEBZ is expected here. )
C
C  ISPLIT  (input) INTEGER array, dimension (N)
C          The splitting points, at which T breaks up into submatrices.
C          The first submatrix consists of rows/columns 1 to
C          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
C          through ISPLIT( 2 ), etc.
C          ( The output array ISPLIT from DSTEBZ is expected here. )
C
C  Z       (output) COMPLEX*16 array, dimension (LDZ, M)
C          The computed eigenvectors.  The eigenvector associated
C          with the eigenvalue W(i) is stored in the i-th column of
C          Z.  Any vector which fails to converge is set to its current
C          iterate after MAXITS iterations.
C
C  LDZ     (input) INTEGER
C          The leading dimension of the array Z.  LDZ >= max(1,N).
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
C
C  IWORK   (workspace) INTEGER array, dimension (N)
C
C  IFAIL   (output) INTEGER array, dimension (M)
C          On normal exit, all elements of IFAIL are zero.
C          If one or more eigenvectors fail to converge after
C          MAXITS iterations, then their indices are stored in
C          array IFAIL.
C
C  INFO    (output) INTEGER
C          = 0: successful exit.
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = +k, then k eigenvectors failed to converge
C               in MAXITS iterations.  Their indices are stored in
C               array IFAIL.
C
C  Internal Parameters
C  ===================
C
C  MAXITS  INTEGER, default = 5
C          The maximum number of iterations performed.
C
C  EXTRA   INTEGER, default = 2
C          The number of iterations performed after norm growth
C          criterion is satisfied, should be at least 1.
C
C-----------------------------------------------------------------------
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=0.0D+0,CONE=1.0D+0)
      DOUBLE PRECISION  ZERO, ONE, TEN, ODM3, ODM1
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TEN=1.0D1,ODM3=1.0D-3,
     *                  ODM1=1.0D-1)
      INTEGER           MAXITS, EXTRA
      PARAMETER         (MAXITS=5,EXTRA=2)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDZ, M, N
C     .. Array Arguments ..
      COMPLEX*16        Z(LDZ,*)
      DOUBLE PRECISION  D(*), E(*), W(*), WORK(*)
      INTEGER           IBLOCK(*), IFAIL(*), ISPLIT(*), IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CTR, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL, SCL,
     *                  SEP, STPCRT, TOL, XJ, XJM
      INTEGER           B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1, INDRV2,
     *                  INDRV3, INDRV4, INDRV5, ITS, J, J1, JBLK, JMAX,
     *                  JR, NBLK, NRMCHK
C     .. Local Arrays ..
      INTEGER           ISEED(4)
C     .. External Functions ..
      DOUBLE PRECISION  DASUM, DNRM2, X02AJF
      INTEGER           IDAMAX, X02BHF
      EXTERNAL          DASUM, DNRM2, X02AJF, IDAMAX, X02BHF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, F01LEF, F04LEF, F06AAZ, G05FAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCMPLX, MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      DO 20 I = 1, M
         IFAIL(I) = 0
   20 CONTINUE
C
      IF (N.LT.0) THEN
         INFO = -1
      ELSE IF (M.LT.0 .OR. M.GT.N) THEN
         INFO = -4
      ELSE IF (LDZ.LT.MAX(1,N)) THEN
         INFO = -9
      ELSE
         DO 40 J = 2, M
            IF (IBLOCK(J).LT.IBLOCK(J-1)) THEN
               INFO = -6
               GO TO 60
            END IF
            IF (IBLOCK(J).EQ.IBLOCK(J-1) .AND. W(J).LT.W(J-1)) THEN
               INFO = -5
               GO TO 60
            END IF
   40    CONTINUE
   60    CONTINUE
      END IF
C
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08JXF/ZSTEIN',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. M.EQ.0) THEN
         RETURN
      ELSE IF (N.EQ.1) THEN
         Z(1,1) = CONE
         RETURN
      END IF
C
C     Get machine constants.
C
      EPS = X02AJF()*X02BHF()
C
C     Initialize seed for random number generator DLARNV.
C
      DO 80 I = 1, 4
         ISEED(I) = 1
   80 CONTINUE
C
C     Initialize pointers.
C
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
C
C     Compute eigenvectors of matrix blocks.
C
      J1 = 1
      DO 360 NBLK = 1, IBLOCK(M)
C
C        Find starting and ending indices of block nblk.
C
         IF (NBLK.EQ.1) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT(NBLK-1) + 1
         END IF
         BN = ISPLIT(NBLK)
         BLKSIZ = BN - B1 + 1
         IF (BLKSIZ.EQ.1) GO TO 120
         GPIND = B1
C
C        Compute reorthogonalization criterion and stopping criterion.
C
         ONENRM = ABS(D(B1)) + ABS(E(B1))
         ONENRM = MAX(ONENRM,ABS(D(BN))+ABS(E(BN-1)))
         DO 100 I = B1 + 1, BN - 1
            ONENRM = MAX(ONENRM,ABS(D(I))+ABS(E(I-1))+ABS(E(I)))
  100    CONTINUE
         ORTOL = ODM3*ONENRM
C
         STPCRT = SQRT(ODM1/BLKSIZ)
C
C        Loop through eigenvalues of block nblk.
C
  120    CONTINUE
         JBLK = 0
         DO 340 J = J1, M
            IF (IBLOCK(J).NE.NBLK) THEN
               J1 = J
               GO TO 360
            END IF
            JBLK = JBLK + 1
            XJ = W(J)
C
C           Skip all the work if the block size is one.
C
            IF (BLKSIZ.EQ.1) THEN
               WORK(INDRV1+1) = ONE
               GO TO 280
            END IF
C
C           If eigenvalues j and j-1 are too close, add a relatively
C           small perturbation.
C
            IF (JBLK.GT.1) THEN
               EPS1 = ABS(EPS*XJ)
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF (SEP.LT.PERTOL) XJ = XJM + PERTOL
            END IF
C
            ITS = 0
            NRMCHK = 0
C
C           Get random starting vector.
C
            CALL G05FAF(-ONE,ONE,BLKSIZ,WORK(INDRV1+1))
C
C           Copy the matrix T so it won't be destroyed in factorization.
C
            CALL DCOPY(BLKSIZ,D(B1),1,WORK(INDRV4+1),1)
            CALL DCOPY(BLKSIZ-1,E(B1),1,WORK(INDRV2+2),1)
            CALL DCOPY(BLKSIZ-1,E(B1),1,WORK(INDRV3+1),1)
C
C           Compute LU factors with partial pivoting  ( PT = LU )
C
            TOL = ZERO
            IINFO = 0
            CALL F01LEF(BLKSIZ,WORK(INDRV4+1),XJ,WORK(INDRV2+2-1),
     *                  WORK(INDRV3+1-1),TOL,WORK(INDRV5+1-2),IWORK,
     *                  IINFO)
C
C           Update iteration count.
C
  140       CONTINUE
            ITS = ITS + 1
            IF (ITS.GT.MAXITS) GO TO 240
C
C           Normalize and scale the righthand side vector Pb.
C
            SCL = BLKSIZ*ONENRM*MAX(EPS,ABS(WORK(INDRV4+BLKSIZ)))
     *            /DASUM(BLKSIZ,WORK(INDRV1+1),1)
            CALL DSCAL(BLKSIZ,SCL,WORK(INDRV1+1),1)
C
C           Solve the system LU = Pb.
C
            IINFO = 0
            CALL F04LEF(-1,BLKSIZ,WORK(INDRV4+1),WORK(INDRV2+2-1),
     *                  WORK(INDRV3+1-1),WORK(INDRV5+1-2),IWORK,
     *                  WORK(INDRV1+1),TOL,IINFO)
C
C           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
C           close enough.
C
            IF (JBLK.EQ.1) GO TO 220
            IF (ABS(XJ-XJM).GT.ORTOL) GPIND = J
            DO 200 I = GPIND, J - 1
               CTR = ZERO
               DO 160 JR = 1, BLKSIZ
                  CTR = CTR + WORK(INDRV1+JR)*DBLE(Z(B1-1+JR,I))
  160          CONTINUE
               DO 180 JR = 1, BLKSIZ
                  WORK(INDRV1+JR) = WORK(INDRV1+JR) -
     *                              CTR*DBLE(Z(B1-1+JR,I))
  180          CONTINUE
  200       CONTINUE
C
C           Check the infinity norm of the iterate.
C
  220       CONTINUE
            JMAX = IDAMAX(BLKSIZ,WORK(INDRV1+1),1)
            NRM = ABS(WORK(INDRV1+JMAX))
C
C           Continue for additional iterations after norm reaches
C           stopping criterion.
C
            IF (NRM.LT.STPCRT) GO TO 140
            NRMCHK = NRMCHK + 1
            IF (NRMCHK.LT.EXTRA+1) GO TO 140
C
            GO TO 260
C
C           If stopping criterion was not satisfied, update info and
C           store eigenvector number in array ifail.
C
  240       CONTINUE
            INFO = INFO + 1
            IFAIL(INFO) = J
C
C           Accept iterate as jth eigenvector.
C
  260       CONTINUE
            SCL = ONE/DNRM2(BLKSIZ,WORK(INDRV1+1),1)
            JMAX = IDAMAX(BLKSIZ,WORK(INDRV1+1),1)
            IF (WORK(INDRV1+JMAX).LT.ZERO) SCL = -SCL
            CALL DSCAL(BLKSIZ,SCL,WORK(INDRV1+1),1)
  280       CONTINUE
            DO 300 I = 1, N
               Z(I,J) = CZERO
  300       CONTINUE
            DO 320 I = 1, BLKSIZ
               Z(B1+I-1,J) = DCMPLX(WORK(INDRV1+I),ZERO)
  320       CONTINUE
C
C           Save the shift to check eigenvalue spacing at next
C           iteration.
C
            XJM = XJ
C
  340    CONTINUE
  360 CONTINUE
C
      RETURN
C
C     End of F08JXF (ZSTEIN)
C
      END
