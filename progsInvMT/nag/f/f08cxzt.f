      SUBROUTINE F08CXZ(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK,
     *                  INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  ZUNMRQ overwrites the general complex M-by-N matrix C with
C
C                  SIDE = 'L'     SIDE = 'R'
C  TRANS = 'N':      Q * C          C * Q
C  TRANS = 'C':      Q**H * C       C * Q**H
C
C  where Q is a complex unitary matrix defined as the product of k
C  elementary reflectors
C
C        Q = H(1)' H(2)' . . . H(k)'
C
C  as returned by ZGERQF. Q is of order M if SIDE = 'L' and of order N
C  if SIDE = 'R'.
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': apply Q or Q**H from the Left;
C          = 'R': apply Q or Q**H from the Right.
C
C  TRANS   (input) CHARACTER*1
C          = 'N':  No transpose, apply Q;
C          = 'C':  Transpose, apply Q**H.
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
C  A       (input) COMPLEX*16 array, dimension
C                               (LDA,M) if SIDE = 'L',
C                               (LDA,N) if SIDE = 'R'
C          The i-th row must contain the vector which defines the
C          elementary reflector H(i), for i = 1,2,...,k, as returned by
C          ZGERQF in the last k rows of its array argument A.
C          A is modified by the routine but restored on exit.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,K).
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGERQF.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the M-by-N matrix C.
C          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDC >= max(1,M).
C
C  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.
C          If SIDE = 'L', LWORK >= max(1,N);
C          if SIDE = 'R', LWORK >= max(1,M).
C          For optimum performance LWORK >= N*NB if SIDE = 'L', and
C          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
C          blocksize.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      INTEGER           NBMAX, LDT
      PARAMETER         (NBMAX=64,LDT=NBMAX+1)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, LDC, LWORK, M, N
      CHARACTER         SIDE, TRANS
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), C(LDC,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, I1, I2, I3, IB, IINFO, IWS, LDWORK, MI, NB,
     *                  NBMIN, NI, NQ, NW
      LOGICAL           LEFT, NOTRAN
      CHARACTER         TRANST
C     .. Local Arrays ..
      COMPLEX*16        T(LDT,NBMAX)
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08ASX, F08ASY, F08CXY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LEFT = SIDE .EQ. 'L' .OR. SIDE .EQ. 'l'
      NOTRAN = TRANS .EQ. 'N' .OR. TRANS .EQ. 'n'
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
      ELSE IF ( .NOT. NOTRAN .AND. .NOT.
     *         (TRANS.EQ.'C' .OR. TRANS.EQ.'c')) THEN
         INFO = -2
      ELSE IF (M.LT.0) THEN
         INFO = -3
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (K.LT.0 .OR. K.GT.NQ) THEN
         INFO = -5
      ELSE IF (LDA.LT.MAX(1,K)) THEN
         INFO = -7
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = -10
      ELSE IF (LWORK.LT.MAX(1,NW)) THEN
         INFO = -12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08CXZ',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size.  NB may be at most NBMAX, where NBMAX
C     is used to define the local array T.
C
      CALL F07ZAZ(1,'F08CXZ',NB,0)
      NB = MIN(NBMAX,NB)
      NBMIN = 2
      LDWORK = NW
      IF (NB.GT.1 .AND. NB.LT.K) THEN
         IWS = NW*NB
         IF (LWORK.LT.IWS) THEN
            NB = LWORK/LDWORK
            CALL F07ZAZ(2,'F08CXZ',NBMIN,0)
            NBMIN = MAX(2,NBMIN)
         END IF
      ELSE
         IWS = NW
      END IF
C
      IF (NB.LT.NBMIN .OR. NB.GE.K) THEN
C
C        Use unblocked code
C
         CALL F08CXY(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,IINFO)
      ELSE
C
C        Use blocked code
C
         IF ((LEFT .AND. .NOT. NOTRAN) .OR. ( .NOT. LEFT .AND. NOTRAN))
     *       THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ((K-1)/NB)*NB + 1
            I2 = 1
            I3 = -NB
         END IF
C
         IF (LEFT) THEN
            NI = N
         ELSE
            MI = M
         END IF
C
         IF (NOTRAN) THEN
            TRANST = 'C'
         ELSE
            TRANST = 'N'
         END IF
C
         DO 20 I = I1, I2, I3
            IB = MIN(NB,K-I+1)
C
C           Form the triangular factor of the block reflector
C           H = H(i+ib-1) . . . H(i+1) H(i)
C
            CALL F08ASX('Backward','Rowwise',NQ-K+I+IB-1,IB,A(I,1),LDA,
     *                  TAU(I),T,LDT)
            IF (LEFT) THEN
C
C              H or H' is applied to C(1:m-k+i+ib-1,1:n)
C
               MI = M - K + I + IB - 1
            ELSE
C
C              H or H' is applied to C(1:m,1:n-k+i+ib-1)
C
               NI = N - K + I + IB - 1
            END IF
C
C           Apply H or H'
C
            CALL F08ASY(SIDE,TRANST,'Backward','Rowwise',MI,NI,IB,A(I,1)
     *                  ,LDA,T,LDT,C,LDC,WORK,LDWORK)
   20    CONTINUE
      END IF
      WORK(1) = IWS
      RETURN
C
C     End of F08CXZ (ZUNMRQ)
C
      END
