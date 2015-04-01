      SUBROUTINE F08GGF(SIDE,UPLO,TRANS,M,N,AP,TAU,C,LDC,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DOPMTR(SIDE,UPLO,TRANS,M,N,AP,TAU,C,LDC,WORK,
     *                  INFO)
C
C  Purpose
C  =======
C
C  DOPMTR overwrites the general real m by n matrix C with
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
C  nq-1 elementary reflectors, as returned by DSPTRD using packed
C  storage:
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
C          Specifies the storage scheme used in the previous call of
C          DSPTRD:
C          = 'U': Upper triangular packed storage;
C          = 'L': Lower triangular packed storage.
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
C  AP      (input) DOUBLE PRECISION array, dimension
C                               (M*(M+1)/2) if SIDE = 'L'
C                               (N*(N+1)/2) if SIDE = 'R'
C          The vectors which define the elementary reflectors, as
C          returned by DSPTRD.
C
C  TAU     (input) DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L'
C                                     or (N-1) if SIDE = 'R'
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by DSPTRD.
C
C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C          On entry, the m by n matrix C.
C          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDC >= max(1,M).
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension
C                                   (N) if SIDE = 'L'
C                                   (M) if SIDE = 'R'
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
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDC, M, N
      CHARACTER         SIDE, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  AP(*), C(LDC,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AII
      INTEGER           I, I1, I2, I3, IC, II, JC, MI, NI, NQ
      LOGICAL           FORWRD, LEFT, NOTRAN, UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08AEW
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LEFT = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
      NOTRAN = (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
C     NQ is the order of Q
C
      IF (LEFT) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF ( .NOT. LEFT .AND. .NOT. (SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
         INFO = -1
      ELSE IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l'))
     *         THEN
         INFO = -2
      ELSE IF ( .NOT. NOTRAN .AND. .NOT.
     *         (TRANS.EQ.'T' .OR. TRANS.EQ.'t')) THEN
         INFO = -3
      ELSE IF (M.LT.0) THEN
         INFO = -4
      ELSE IF (N.LT.0) THEN
         INFO = -5
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = -9
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08GGF/DOPMTR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Q was determined by a call to DSPTRD with UPLO = 'U'
C
         FORWRD = (LEFT .AND. NOTRAN) .OR.
     *            ( .NOT. LEFT .AND. .NOT. NOTRAN)
C
         IF (FORWRD) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*(NQ+1)/2 - 1
         END IF
C
         IF (LEFT) THEN
            NI = N
         ELSE
            MI = M
         END IF
C
         DO 20 I = I1, I2, I3
            IF (LEFT) THEN
C
C              H(i) is applied to C(1:i,1:n)
C
               MI = I
            ELSE
C
C              H(i) is applied to C(1:m,1:i)
C
               NI = I
            END IF
C
C           Apply H(i)
C
            AII = AP(II)
            AP(II) = ONE
            CALL F08AEW(SIDE,MI,NI,AP(II-I+1),1,TAU(I),C,LDC,WORK)
            AP(II) = AII
C
            IF (FORWRD) THEN
               II = II + I + 2
            ELSE
               II = II - I - 1
            END IF
   20    CONTINUE
      ELSE
C
C        Q was determined by a call to DSPTRD with UPLO = 'L'.
C
         FORWRD = (LEFT .AND. .NOT. NOTRAN)
     *            .OR. ( .NOT. LEFT .AND. NOTRAN)
C
         IF (FORWRD) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*(NQ+1)/2 - 1
         END IF
C
         IF (LEFT) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
C
         DO 40 I = I1, I2, I3
            AII = AP(II)
            AP(II) = ONE
            IF (LEFT) THEN
C
C              H(i) is applied to C(i+1:m,1:n)
C
               MI = M - I
               IC = I + 1
            ELSE
C
C              H(i) is applied to C(1:m,i+1:n)
C
               NI = N - I
               JC = I + 1
            END IF
C
C           Apply H(i)
C
            CALL F08AEW(SIDE,MI,NI,AP(II),1,TAU(I),C(IC,JC),LDC,WORK)
            AP(II) = AII
C
            IF (FORWRD) THEN
               II = II + NQ - I + 1
            ELSE
               II = II - NQ + I - 2
            END IF
   40    CONTINUE
      END IF
      RETURN
C
C     End of F08GGF (DOPMTR)
C
      END
