      SUBROUTINE F08KGF(VECT,SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,
     *                  LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DORMBR(VECT,SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,
     *                  WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  If VECT = 'Q', DORMBR overwrites the general real m by n matrix C
C  with
C
C        Q * C  if SIDE = 'L' and TRANS = 'N', or
C
C        Q'* C  if SIDE = 'L' and TRANS = 'T', or
C
C        C * Q  if SIDE = 'R' and TRANS = 'N', or
C
C        C * Q' if SIDE = 'R' and TRANS = 'T'.
C
C  If VECT = 'P', DORMBR overwrites the general real m by n matrix C
C  with
C
C        P * C  if SIDE = 'L' and TRANS = 'N', or
C
C        P'* C  if SIDE = 'L' and TRANS = 'T', or
C
C        C * P  if SIDE = 'R' and TRANS = 'N', or
C
C        C * P' if SIDE = 'R' and TRANS = 'T'.
C
C  Here Q and P' are the orthogonal matrices determined by DGEBRD when
C  reducing a real matrix A to bidiagonal form: A = Q * B * P'. Q and P'
C  are defined as products of elementary reflectors H(i) and G(i)
C  respectively.
C
C  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
C  order of the orthogonal matrix Q or P' that is applied.
C
C  If VECT = 'Q', A is assumed to have been an nq-by-k matrix:
C  if nq >= k, Q = H(1) H(2) . . . H(k);
C  if nq < k, Q = H(1) H(2) . . . H(nq-1).
C
C  If VECT = 'P', A is assumed to have been a k-by-nq matrix:
C  if k < nq, P = G(1) G(2) . . . G(k);
C  if k >= nq, P = G(1) G(2) . . . G(nq-1).
C
C  Arguments
C  =========
C
C  VECT    (input) CHARACTER*1
C          = 'Q': apply Q or Q'
C          = 'P': apply P or P'
C
C  SIDE    (input) CHARACTER*1
C          = 'L': apply Q, Q', P or P' from the Left
C          = 'R': apply Q, Q', P or P' from the Right
C
C  TRANS   (input) CHARACTER*1
C          = 'N': apply Q  or P  (No transpose)
C          = 'T': apply Q' or P' (Transpose)
C
C  M       (input) INTEGER
C          The number of rows of the matrix C. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C. N >= 0.
C
C  K       (input) INTEGER
C          If VECT = 'Q', the number of columns in the original
C          matrix reduced by DGEBRD.
C          If VECT = 'P', the number of rows in the original
C          matrix reduced by DGEBRD.
C          K >= 0.
C
C  A       (input) DOUBLE PRECISION array, dimension
C                                (LDA,min(nq,K)) if VECT = 'Q'
C                                (LDA,nq)        if VECT = 'P'
C          The vectors which define the elementary reflectors H(i) and
C          G(i), whose products determine the matrices Q and P, as
C          returned by DGEBRD.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.
C          If VECT = 'Q', LDA >= max(1,nq);
C          if VECT = 'P', LDA >= max(1,min(nq,K)).
C
C  TAU     (input) DOUBLE PRECISION array, dimension (min(nq,K))
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i) or G(i) which determines Q or P, as returned
C          by DGEBRD in the array argument TAUQ or TAUP.
C
C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C          On entry, the m by n matrix C.
C          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q
C          or P*C or P'*C or C*P or C*P'.
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
      INTEGER           INFO, K, LDA, LDC, LWORK, M, N
      CHARACTER         SIDE, TRANS, VECT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(LDC,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I1, I2, IINFO, MI, NI, NQ, NW
      LOGICAL           APPLYQ, LEFT, NOTRAN
      CHARACTER         TRANST
C     .. External Subroutines ..
      EXTERNAL          DORMLQ, DORMQR, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      APPLYQ = (VECT.EQ.'Q' .OR. VECT.EQ.'q')
      LEFT = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
      NOTRAN = (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
C
C     NQ is the order of Q or P and NW is the minimum dimension of WORK
C
      IF (LEFT) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF ( .NOT. APPLYQ .AND. .NOT. (VECT.EQ.'P' .OR. VECT.EQ.'p')) THEN
         INFO = -1
      ELSE IF ( .NOT. LEFT .AND. .NOT. (SIDE.EQ.'R' .OR. SIDE.EQ.'r'))
     *         THEN
         INFO = -2
      ELSE IF ( .NOT. NOTRAN .AND. .NOT.
     *         (TRANS.EQ.'T' .OR. TRANS.EQ.'t')) THEN
         INFO = -3
      ELSE IF (M.LT.0) THEN
         INFO = -4
      ELSE IF (N.LT.0) THEN
         INFO = -5
      ELSE IF (K.LT.0) THEN
         INFO = -6
      ELSE IF ((APPLYQ .AND. LDA.LT.MAX(1,NQ))
     *         .OR. ( .NOT. APPLYQ .AND. LDA.LT.MAX(1,MIN(NQ,K)))) THEN
         INFO = -8
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = -11
      ELSE IF (LWORK.LT.MAX(1,NW)) THEN
         INFO = -13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08KGF/DORMBR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
      IF (APPLYQ) THEN
C
C        Apply Q
C
         IF (NQ.GE.K) THEN
C
C           Q was determined by a call to DGEBRD with nq >= k
C
            CALL DORMQR(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK,
     *                  IINFO)
         ELSE IF (NQ.GT.1) THEN
C
C           Q was determined by a call to DGEBRD with nq < k
C
            IF (LEFT) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL DORMQR(SIDE,TRANS,MI,NI,NQ-1,A(2,1),LDA,TAU,C(I1,I2),
     *                  LDC,WORK,LWORK,IINFO)
         ELSE
            WORK(1) = 1
         END IF
      ELSE
C
C        Apply P
C
         IF (NOTRAN) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
         IF (NQ.GT.K) THEN
C
C           P was determined by a call to DGEBRD with nq > k
C
            CALL DORMLQ(SIDE,TRANST,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK,
     *                  IINFO)
         ELSE IF (NQ.GT.1) THEN
C
C           P was determined by a call to DGEBRD with nq <= k
C
            IF (LEFT) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL DORMLQ(SIDE,TRANST,MI,NI,NQ-1,A(1,2),LDA,TAU,C(I1,I2),
     *                  LDC,WORK,LWORK,IINFO)
         ELSE
            WORK(1) = 1
         END IF
      END IF
      RETURN
C
C     End of F08KGF (DORMBR)
C
      END
