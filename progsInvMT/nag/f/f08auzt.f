      SUBROUTINE F08AUZ(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZUNM2R(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,
C    *                  INFO)
C
C  Purpose
C  =======
C
C  ZUNM2R overwrites the general complex m-by-n matrix C with
C
C        Q * C  if SIDE = 'L' and TRANS = 'N', or
C
C        Q'* C  if SIDE = 'L' and TRANS = 'C', or
C
C        C * Q  if SIDE = 'R' and TRANS = 'N', or
C
C        C * Q' if SIDE = 'R' and TRANS = 'C',
C
C  where Q is a complex unitary matrix defined as the product of k
C  elementary reflectors
C
C        Q = H(1) H(2) . . . H(k)
C
C  as returned by ZGEQRF. Q is of order m if SIDE = 'L' and of order n
C  if SIDE = 'R'.
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
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines
C          the matrix Q.
C          If SIDE = 'L', M >= K >= 0;
C          if SIDE = 'R', N >= K >= 0.
C
C  A       (input) COMPLEX*16 array, dimension (LDA,K)
C          The i-th column must contain the vector which defines the
C          elementary reflector H(i), for i = 1,2,...,k, as returned by
C          ZGEQRF in the first k columns of its array argument A.
C          A is modified by the routine but restored on exit.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.
C          If SIDE = 'L', LDA >= max(1,M);
C          if SIDE = 'R', LDA >= max(1,N).
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEQRF.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the m-by-n matrix C.
C          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDC >= max(1,M).
C
C  WORK    (workspace) COMPLEX*16 array, dimension
C                                   (N) if SIDE = 'L',
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
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, LDC, M, N
      CHARACTER         SIDE, TRANS
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), C(LDC,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        AII, TAUI
      INTEGER           I, I1, I2, I3, IC, JC, MI, NI, NQ
      LOGICAL           LEFT, NOTRAN
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASW
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LEFT = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
      NOTRAN = (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
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
      ELSE IF ( .NOT. NOTRAN .AND. .NOT.
     *         (TRANS.EQ.'C' .OR. TRANS.EQ.'c')) THEN
         INFO = -2
      ELSE IF (M.LT.0) THEN
         INFO = -3
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (K.LT.0 .OR. K.GT.NQ) THEN
         INFO = -5
      ELSE IF (LDA.LT.MAX(1,NQ)) THEN
         INFO = -7
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = -10
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08AUZ/ZUNM2R',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0) RETURN
C
      IF ((LEFT .AND. .NOT. NOTRAN .OR. .NOT. LEFT .AND. NOTRAN)) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
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
      DO 20 I = I1, I2, I3
         IF (LEFT) THEN
C
C           H(i) or H(i)' is applied to C(i:m,1:n)
C
            MI = M - I + 1
            IC = I
         ELSE
C
C           H(i) or H(i)' is applied to C(1:m,i:n)
C
            NI = N - I + 1
            JC = I
         END IF
C
C        Apply H(i) or H(i)'
C
         IF (NOTRAN) THEN
            TAUI = TAU(I)
         ELSE
            TAUI = DCONJG(TAU(I))
         END IF
         AII = A(I,I)
         A(I,I) = ONE
         CALL F08ASW(SIDE,MI,NI,A(I,I),1,TAUI,C(IC,JC),LDC,WORK)
         A(I,I) = AII
   20 CONTINUE
      RETURN
C
C     End of F08AUZ (ZUNM2R)
C
      END
