      SUBROUTINE F08CVZ(M,N,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  ZGERQF computes an RQ factorization of a complex M-by-N matrix A:
C  A = R * Q.
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
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the M-by-N matrix A.
C          On exit,
C          if m <= n, the upper triangle of the subarray
C          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
C          if m >= n, the elements on and above the (m-n)-th subdiagonal
C          contain the M-by-N upper trapezoidal matrix R;
C          the remaining elements, with the array TAU, represent the
C          unitary matrix Q as a product of min(m,n) elementary
C          reflectors (see Further Details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.  LWORK >= max(1,M).
C          For optimum performance LWORK >= M*NB, where NB is
C          the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value
C
C  Further Details
C  ===============
C
C  The matrix Q is represented as a product of elementary reflectors
C
C     Q = H(1)' H(2)' . . . H(k)', where k = min(m,n).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on
C  exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i).
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IB, IINFO, IWS, K, KI, KK, LDWORK, MU, NB,
     *                  NBMIN, NU, NX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08ASX, F08ASY, F08CVY
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
      ELSE IF (LWORK.LT.MAX(1,M)) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08CVZ',-INFO)
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
      CALL F07ZAZ(1,'F08CVZ',NB,0)
      NBMIN = 2
      NX = 1
      IWS = M
      IF (NB.GT.1 .AND. NB.LT.K) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         CALL F07ZAZ(3,'F08CVZ',NX,0)
         NX = MAX(0,NX)
         IF (NX.LT.K) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = M
            IWS = LDWORK*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK/LDWORK
               CALL F07ZAZ(2,'F08CVZ',NBMIN,0)
               NBMIN = MAX(2,NBMIN)
            END IF
         END IF
      END IF
C
      IF (NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K) THEN
C
C        Use blocked code initially.
C        The last kk rows are handled by the block method.
C
         KI = ((K-NX-1)/NB)*NB
         KK = MIN(K,KI+NB)
C
         DO 20 I = K - KK + KI + 1, K - KK + 1, -NB
            IB = MIN(K-I+1,NB)
C
C           Compute the RQ factorization of the current block
C           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)
C
            CALL F08CVY(IB,N-K+I+IB-1,A(M-K+I,1),LDA,TAU(I),WORK,IINFO)
            IF (M-K+I.GT.1) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i+ib-1) . . . H(i+1) H(i)
C
               CALL F08ASX('Backward','Rowwise',N-K+I+IB-1,IB,A(M-K+I,1)
     *                     ,LDA,TAU(I),WORK,LDWORK)
C
C              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
C
               CALL F08ASY('Right','No transpose','Backward','Rowwise',
     *                     M-K+I-1,N-K+I+IB-1,IB,A(M-K+I,1),LDA,WORK,
     *                     LDWORK,A,LDA,WORK(IB+1),LDWORK)
            END IF
   20    CONTINUE
         MU = M - K + I + NB - 1
         NU = N - K + I + NB - 1
      ELSE
         MU = M
         NU = N
      END IF
C
C     Use unblocked code to factor the last or only block
C
      IF (MU.GT.0 .AND. NU.GT.0) CALL F08CVY(MU,NU,A,LDA,TAU,WORK,IINFO)
C
      WORK(1) = IWS
      RETURN
C
C     End of F08CVZ (ZGERQF)
C
      END
