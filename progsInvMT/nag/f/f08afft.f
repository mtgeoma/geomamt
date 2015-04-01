      SUBROUTINE F08AFF(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DORGQR(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  DORGQR generates an m by n real matrix Q with orthonormal columns,
C  which is defined as the first n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(1) H(2) . . . H(k)
C
C  as returned by DGEQRF.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q. M >= N >= 0.
C
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines the
C          matrix Q. N >= K >= 0.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the i-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by DGEQRF in the first k columns of its array
C          argument A.
C          On exit, the m by n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) DOUBLE PRECISION array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by DGEQRF.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,N).
C          For optimum performance LWORK should be at least N*NB, where
C          NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument has an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, LWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, NB,
     *                  NBMIN, NX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08AEX, F08AEY, F08AFZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0 .OR. N.GT.M) THEN
         INFO = -2
      ELSE IF (K.LT.0 .OR. K.GT.N) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      ELSE IF (LWORK.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08AFF/DORGQR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size.
C
      CALL F07ZAZ(1,'F08AFF',NB,0)
      NBMIN = 2
      NX = 0
      IWS = N
      IF (NB.GT.1 .AND. NB.LT.K) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         CALL F07ZAZ(3,'F08AFF',NX,0)
         NX = MAX(0,NX)
         IF (NX.LT.K) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK/LDWORK
               CALL F07ZAZ(2,'F08AFF',NBMIN,0)
               NBMIN = MAX(2,NBMIN)
            END IF
         END IF
      END IF
C
      IF (NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K) THEN
C
C        Use blocked code after the last block.
C        The first kk columns are handled by the block method.
C
         KI = ((K-NX-1)/NB)*NB
         KK = MIN(K,KI+NB)
C
C        Set A(1:kk,kk+1:n) to zero.
C
         DO 40 J = KK + 1, N
            DO 20 I = 1, KK
               A(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the last or only block.
C
      IF (KK.LT.N) CALL F08AFZ(M-KK,N-KK,K-KK,A(KK+1,KK+1),LDA,TAU(KK+1)
     *                         ,WORK,IINFO)
C
      IF (KK.GT.0) THEN
C
C        Use blocked code
C
         DO 100 I = KI + 1, 1, -NB
            IB = MIN(NB,K-I+1)
            IF (I+IB.LE.N) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i) H(i+1) . . . H(i+ib-1)
C
               CALL F08AEX('Forward','Columnwise',M-I+1,IB,A(I,I),LDA,
     *                     TAU(I),WORK,LDWORK)
C
C              Apply H to A(i:m,i+ib:n) from the left
C
               CALL F08AEY('Left','No transpose','Forward','Columnwise',
     *                     M-I+1,N-I-IB+1,IB,A(I,I),LDA,WORK,LDWORK,
     *                     A(I,I+IB),LDA,WORK(IB+1),LDWORK)
            END IF
C
C           Apply H to rows i:m of current block
C
            CALL F08AFZ(M-I+1,IB,IB,A(I,I),LDA,TAU(I),WORK,IINFO)
C
C           Set rows 1:i-1 of current block to zero
C
            DO 80 J = I, I + IB - 1
               DO 60 L = 1, I - 1
                  A(L,J) = ZERO
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
      END IF
C
      WORK(1) = IWS
      RETURN
C
C     End of F08AFF (DORGQR)
C
      END
