      SUBROUTINE F08BEF(M,N,A,LDA,JPVT,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DGEQPF(M,N,A,LDA,JPVT,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  DGEQPF computes a QR factorization with column pivoting of a
C  real m by n matrix A: A*P = Q*R
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A. N >= 0
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the m by n matrix A.
C          On exit, the upper triangle of the array contains the
C          min(m,n) by n upper triangular matrix R; the elements
C          below the diagonal, together with the array TAU,
C          represent the orthogonal matrix Q as a product of
C          min(m,n) elementary reflectors.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,M).
C
C  JPVT    (input/output) INTEGER array, dimension (N)
C          on entry: If JPVT(I) <> 0, column I of A is permuted
C          to the front of AP (a leading column)
C          IF JPVT(I) == 0, column I of A is a free column.
C          on exit: If JPVT(I) = K, then the Ith column of AP
C          was the Kth column of A.
C
C  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
C          Stores further details of
C          the orthogonal matrix Q (see A).
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
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
C     Q = H(1) H(2) . . . H(n)
C
C  Each H(i) has the form
C
C     H = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with
C  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
C
C  The matrix P is represented in jpvt as follows: If
C     jpvt(j) = i
C  then the jth column of P is the ith canonical unit vector.
C
C  -- LAPACK test routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), TAU(*), WORK(*)
      INTEGER           JPVT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AII, TEMP, TEMP2
      INTEGER           I, ITEMP, J, MA, MN, PVT
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DSWAP, F06AAZ, F08AEV, F08AEW, F08AEZ, F08AGZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
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
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08BEF/DGEQPF',-INFO)
         RETURN
      END IF
C
      MN = MIN(M,N)
C
C     Move initial columns up front
C
      ITEMP = 1
      DO 20 I = 1, N
         IF (JPVT(I).NE.0) THEN
            IF (I.NE.ITEMP) THEN
               CALL DSWAP(M,A(1,I),1,A(1,ITEMP),1)
               JPVT(I) = JPVT(ITEMP)
               JPVT(ITEMP) = I
            ELSE
               JPVT(I) = I
            END IF
            ITEMP = ITEMP + 1
         ELSE
            JPVT(I) = I
         END IF
   20 CONTINUE
      ITEMP = ITEMP - 1
C
C     Compute the QR factorization and update remaining columns
C
      IF (ITEMP.GT.0) THEN
         MA = MIN(ITEMP,M)
         CALL F08AEZ(M,MA,A,LDA,TAU,WORK,INFO)
         IF (MA.LT.N) THEN
            CALL F08AGZ('Left','Transpose',M,N-MA,MA,A,LDA,TAU,A(1,MA+1)
     *                  ,LDA,WORK,INFO)
         END IF
      END IF
C
      IF (ITEMP.LT.MN) THEN
C
C        Initialize partial column norms. The first n entries of
C        work store the exact column norms.
C
         DO 40 I = ITEMP + 1, N
            WORK(I) = DNRM2(M-ITEMP,A(ITEMP+1,I),1)
            WORK(N+I) = WORK(I)
   40    CONTINUE
C
C        Compute factorization
C
         DO 80 I = ITEMP + 1, MN
C
C           Determine ith pivot column and swap if necessary
C
            PVT = (I-1) + IDAMAX(N-I+1,WORK(I),1)
C
            IF (PVT.NE.I) THEN
               CALL DSWAP(M,A(1,PVT),1,A(1,I),1)
               ITEMP = JPVT(PVT)
               JPVT(PVT) = JPVT(I)
               JPVT(I) = ITEMP
               WORK(PVT) = WORK(I)
               WORK(N+PVT) = WORK(N+I)
            END IF
C
C           Generate elementary reflector H(i)
C
            IF (I.LT.M) THEN
               CALL F08AEV(M-I+1,A(I,I),A(I+1,I),1,TAU(I))
            ELSE
               CALL F08AEV(1,A(M,M),A(M,M),1,TAU(M))
            END IF
C
            IF (I.LT.N) THEN
C
C              Apply H(i) to A(i:m,i+1:n) from the left
C
               AII = A(I,I)
               A(I,I) = ONE
               CALL F08AEW('LEFT',M-I+1,N-I,A(I,I),1,TAU(I),A(I,I+1),
     *                     LDA,WORK(2*N+1))
               A(I,I) = AII
            END IF
C
C           Update partial column norms
C
            DO 60 J = I + 1, N
               IF (WORK(J).NE.ZERO) THEN
                  TEMP = ONE - (ABS(A(I,J))/WORK(J))**2
                  TEMP = MAX(TEMP,ZERO)
                  TEMP2 = ONE + 0.05D0*TEMP*(WORK(J)/WORK(N+J))**2
                  IF (TEMP2.EQ.ONE) THEN
                     IF (M-I.GT.0) THEN
                        WORK(J) = DNRM2(M-I,A(I+1,J),1)
                        WORK(N+J) = WORK(J)
                     ELSE
                        WORK(J) = ZERO
                        WORK(N+J) = ZERO
                     END IF
                  ELSE
                     WORK(J) = WORK(J)*SQRT(TEMP)
                  END IF
               END IF
   60       CONTINUE
C
   80    CONTINUE
      END IF
      RETURN
C
C     End of F08BEF (DGEQPF)
C
      END
