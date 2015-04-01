      SUBROUTINE F08FGF(SIDE,UPLO,TRANS,M,N,A,LDA,TAU,C,LDC,WORK,LWORK,
     *                  INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DORMTR(SIDE,UPLO,TRANS,M,N,A,LDA,TAU,C,LDC,WORK,
     *                  LWORK,INFO)
C
C  Purpose
C  =======
C
C  DORMTR overwrites the general real m by n matrix C with
C
C        Q * C  if SIDE = 'L' and TRANS = 'N', or
C
C        Q'* C  if SIDE = 'L' and TRANS = 'T', or
C
C        C * Q  if SIDE = 'R' and TRANS = 'N', or
C
C        C * Q' if SIDE = 'R' and TRANS = 'T',
C
C  where Q is a real orthogonal matrix of order nq, with nq = m if
C  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
C  nq-1 elementary reflectors, as returned by DSYTRD:
C
C  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
C
C  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': apply Q or Q' from the Left
C          = 'R': apply Q or Q' from the Right
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangle of the array A
C          holds details of the elementary reflectors, as returned by
C          DSYTRD:
C          = 'U': Upper triangle;
C          = 'L': Lower triangle.
C
C  TRANS   (input) CHARACTER*1
C          = 'N': apply Q  (No transpose)
C          = 'T': apply Q' (Transpose)
C
C  M       (input) INTEGER
C          The number of rows of the matrix C. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C. N >= 0.
C
C  A       (input) DOUBLE PRECISION array, dimension
C                               (LDA,M) if SIDE = 'L'
C                               (LDA,N) if SIDE = 'R'
C          The vectors which define the elementary reflectors, as
C          returned by DSYTRD.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.
C          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
C
C  TAU     (input) DOUBLE PRECISION array, dimension
C                               (M-1) if SIDE = 'L'
C                               (N-1) if SIDE = 'R'
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by DSYTRD.
C
C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C          On entry, the m-by-n matrix C.
C          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDC >= max(1,M).
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.
C          If SIDE = 'L', LWORK >= max(1,N);
C          if SIDE = 'R', LWORK >= max(1,M).
C          For optimum performance LWORK should be at least N*NB
C          if SIDE = 'L' and at least M*NB if SIDE = 'R', where NB is
C          the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDC, LWORK, M, N
      CHARACTER         SIDE, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(LDC,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I1, I2, IINFO, MI, NI, NQ, NW
      LOGICAL           LEFT, UPPER
C     .. External Subroutines ..
      EXTERNAL          DORMQR, F06AAZ, F08CGY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LEFT = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
C     NQ is the order of Q and NW is the minimum dimension of WORK
C
      IF (LEFT) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF ( .NOT. LEFT .AND. .NOT. (SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
         INFO = -1
      ELSE IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l'))
     *         THEN
         INFO = -2
      ELSE IF ( .NOT. (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
     *         .AND. .NOT. (TRANS.EQ.'T' .OR. TRANS.EQ.'t')) THEN
         INFO = -3
      ELSE IF (M.LT.0) THEN
         INFO = -4
      ELSE IF (N.LT.0) THEN
         INFO = -5
      ELSE IF (LDA.LT.MAX(1,NQ)) THEN
         INFO = -7
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = -10
      ELSE IF (LWORK.LT.MAX(1,NW)) THEN
         INFO = -12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08FGF/DORMTR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1) THEN
         WORK(1) = 1
         RETURN
      END IF
C
      IF (LEFT) THEN
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      END IF
C
      IF (UPPER) THEN
C
C        Q was determined by a call to DSYTRD with UPLO = 'U'
C
         CALL F08CGY(SIDE,TRANS,MI,NI,NQ-1,A(1,2),LDA,TAU,C,LDC,WORK,
     *               LWORK,IINFO)
      ELSE
C
C        Q was determined by a call to DSYTRD with UPLO = 'L'
C
         IF (LEFT) THEN
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         END IF
         CALL DORMQR(SIDE,TRANS,MI,NI,NQ-1,A(2,1),LDA,TAU,C(I1,I2),LDC,
     *               WORK,LWORK,IINFO)
      END IF
      RETURN
C
C     End of F08FGF (DORMTR)
C
      END
