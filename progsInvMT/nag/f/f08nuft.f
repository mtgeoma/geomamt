      SUBROUTINE F08NUF(SIDE,TRANS,M,N,ILO,IHI,A,LDA,TAU,C,LDC,WORK,
     *                  LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZUNMHR(SIDE,TRANS,M,N,ILO,IHI,A,LDA,TAU,C,LDC,
     *                  WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNMHR overwrites the general complex m-by-n matrix C with
C
C        Q * C  if SIDE = 'L' and TRANS = 'N', or
C
C        Q'* C  if SIDE = 'L' and TRANS = 'C', or
C
C        C * Q  if SIDE = 'R' and TRANS = 'N', or
C
C        C * Q' if SIDE = 'R' and TRANS = 'C',
C
C  where Q is a complex unitary matrix of order nq, with nq = m if
C  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
C  ihi-ilo elementary reflectors, as returned by ZGEHRD:
C
C  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': apply Q or Q' from the Left
C          = 'R': apply Q or Q' from the Right
C
C  TRANS   (input) CHARACTER*1
C          = 'N': apply Q  (No transpose)
C          = 'C': apply Q' (Conjugate transpose)
C
C  M       (input) INTEGER
C          The number of rows of the matrix C. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C. N >= 0.
C
C  ILO     (input) INTEGER
C  IHI     (input) INTEGER
C          ILO and IHI must have the same values as in the previous call
C          of ZGEHRD. Q is equal to the unit matrix except in the
C          submatrix Q(ilo+1:ihi,ilo+1:ihi).
C          If SIDE = 'L', 1 <= ILO <= IHI <= max(1,M);
C          if SIDE = 'R', 1 <= ILO <= IHI <= max(1,N);
C
C  A       (input) COMPLEX*16 array, dimension
C                               (LDA,M) if SIDE = 'L'
C                               (LDA,N) if SIDE = 'R'
C          The vectors which define the elementary reflectors, as
C          returned by ZGEHRD.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.
C          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
C
C  TAU     (input) COMPLEX*16 array, dimension
C                               (M-1) if SIDE = 'L'
C                               (N-1) if SIDE = 'R'
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEHRD.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the m-by-n matrix C.
C          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDC >= max(1,M).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
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
      INTEGER           IHI, ILO, INFO, LDA, LDC, LWORK, M, N
      CHARACTER         SIDE, TRANS
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), C(LDC,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I1, I2, IINFO, MI, NH, NI, NQ, NW
      LOGICAL           LEFT
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZUNMQR
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LEFT = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
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
      ELSE IF ( .NOT. (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
     *         .AND. .NOT. (TRANS.EQ.'C' .OR. TRANS.EQ.'c')) THEN
         INFO = -2
      ELSE IF (M.LT.0) THEN
         INFO = -3
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (ILO.LT.1 .OR. ILO.GT.MAX(1,NQ)) THEN
         INFO = -5
      ELSE IF (IHI.LT.MIN(ILO,NQ) .OR. IHI.GT.NQ) THEN
         INFO = -6
      ELSE IF (LDA.LT.MAX(1,NQ)) THEN
         INFO = -8
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = -11
      ELSE IF (LWORK.LT.MAX(1,NW)) THEN
         INFO = -13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08NUF/ZUNMHR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      NH = IHI - ILO
      IF (M.EQ.0 .OR. N.EQ.0 .OR. NH.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
      IF (LEFT) THEN
         MI = NH
         NI = N
         I1 = ILO + 1
         I2 = 1
      ELSE
         MI = M
         NI = NH
         I1 = 1
         I2 = ILO + 1
      END IF
C
      CALL ZUNMQR(SIDE,TRANS,MI,NI,NH,A(ILO+1,ILO),LDA,TAU(ILO),C(I1,I2)
     *            ,LDC,WORK,LWORK,IINFO)
      RETURN
C
C     End of F08NUF (ZUNMHR)
C
      END
