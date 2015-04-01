      SUBROUTINE F08CFY(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DORGQL(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  DORGQL generates an m by n real matrix Q with orthonormal columns,
C  which is defined as the last n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(k) . . . H(2) H(1)
C
C  as returned by DGEQLF.
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
C          On entry, the (n-k+i)-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by DGEQLF in the last k columns of its array
C          argument A.
C          On exit, the m by n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) DOUBLE PRECISION array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by DGEQLF.
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
      INTEGER           I, IB, IINFO, IWS, J, KK, L, LDWORK, NB, NBMIN,
     *                  NX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08AEX, F08AEY, F08CFZ
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
         CALL F06AAZ('F08CFY/DORGQL',-INFO)
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
      CALL F07ZAZ(1,'F08CFY',NB,0)
      NBMIN = 2
      NX = 0
      IWS = N
      IF (NB.GT.1 .AND. NB.LT.K) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         CALL F07ZAZ(3,'F08CFY',NX,0)
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
               CALL F07ZAZ(2,'F08CFY',NBMIN,0)
               NBMIN = MAX(2,NBMIN)
            END IF
         END IF
      END IF
C
      IF (NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K) THEN
C
C        Use blocked code after the first block.
C        The last kk columns are handled by the block method.
C
         KK = MIN(K,((K-NX+NB-1)/NB)*NB)
C
C        Set A(m-kk+1:m,1:n-kk) to zero.
C
         DO 40 J = 1, N - KK
            DO 20 I = M - KK + 1, M
               A(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the first or only block.
C
      CALL F08CFZ(M-KK,N-KK,K-KK,A,LDA,TAU,WORK,IINFO)
C
      IF (KK.GT.0) THEN
C
C        Use blocked code
C
         DO 100 I = K - KK + 1, K, NB
            IB = MIN(NB,K-I+1)
            IF (N-K+I.GT.1) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i+ib-1) . . . H(i+1) H(i)
C
               CALL F08AEX('Backward','Columnwise',M-K+I+IB-1,IB,
     *                     A(1,N-K+I),LDA,TAU(I),WORK,LDWORK)
C
C              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
C
               CALL F08AEY('Left','No transpose','Backward',
     *                     'Columnwise',M-K+I+IB-1,N-K+I-1,IB,A(1,N-K+I)
     *                     ,LDA,WORK,LDWORK,A,LDA,WORK(IB+1),LDWORK)
            END IF
C
C           Apply H to rows 1:m-k+i+ib-1 of current block
C
            CALL F08CFZ(M-K+I+IB-1,IB,IB,A(1,N-K+I),LDA,TAU(I),WORK,
     *                  IINFO)
C
C           Set rows m-k+i+ib:m of current block to zero
C
            DO 80 J = N - K + I, N - K + I + IB - 1
               DO 60 L = M - K + I + IB, M
                  A(L,J) = ZERO
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
      END IF
C
      WORK(1) = IWS
      RETURN
C
C     End of F08CFY (DORGQL)
C
      END
