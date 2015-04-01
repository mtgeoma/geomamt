      SUBROUTINE F08AEF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  DGEQRF computes a QR factorization of a real m by n matrix A:
C  A = Q * R.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the m by n matrix A.
C          On exit, the elements on and above the diagonal of the array
C          contain the min(m,n) by n upper trapezoidal matrix R (R is
C          upper triangular if m >= n); the elements below the diagonal,
C          with the array TAU, represent the orthogonal matrix Q as a
C          product of elementary reflectors (see Further Details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.  LWORK >= max(1,N).
C          For optimum performance LWORK should be at least N*NB,
C          where NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  Further Details
C  ===============
C
C  The matrix Q is represented as a product of elementary reflectors
C
C     Q = H(1) H(2) . . . H(k), where k = min(m,n).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with
C  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
C  and tau in TAU(i).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IB, IINFO, IWS, K, LDWORK, NB, NBMIN, NX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08AEX, F08AEY, F08AEZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.MAX(1,N)) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08AEF/DGEQRF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      K = MIN(M,N)
      IF (K.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size.
C
      CALL F07ZAZ(1,'F08AEF',NB,0)
      NBMIN = 2
      NX = 0
      IWS = N
      IF (NB.GT.1 .AND. NB.LT.K) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         CALL F07ZAZ(3,'F08AEF',NX,0)
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
               CALL F07ZAZ(2,'F08AEF',NBMIN,0)
               NBMIN = MAX(2,NBMIN)
            END IF
         END IF
      END IF
C
      IF (NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K) THEN
C
C        Use blocked code initially
C
         DO 20 I = 1, K - NX, NB
            IB = MIN(K-I+1,NB)
C
C           Compute the QR factorization of the current block
C           A(i:m,i:i+ib-1)
C
            CALL F08AEZ(M-I+1,IB,A(I,I),LDA,TAU(I),WORK,IINFO)
            IF (I+IB.LE.N) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i) H(i+1) . . . H(i+ib-1)
C
               CALL F08AEX('Forward','Columnwise',M-I+1,IB,A(I,I),LDA,
     *                     TAU(I),WORK,LDWORK)
C
C              Apply H' to A(i:m,i+ib:n) from the left
C
               CALL F08AEY('Left','Transpose','Forward','Columnwise',
     *                     M-I+1,N-I-IB+1,IB,A(I,I),LDA,WORK,LDWORK,
     *                     A(I,I+IB),LDA,WORK(IB+1),LDWORK)
            END IF
   20    CONTINUE
      ELSE
         I = 1
      END IF
C
C     Use unblocked code to factor the last or only block.
C
      IF (I.LE.K) CALL F08AEZ(M-I+1,N-I+1,A(I,I),LDA,TAU(I),WORK,IINFO)
C
      WORK(1) = IWS
      RETURN
C
C     End of F08AEF (DGEQRF)
C
      END
