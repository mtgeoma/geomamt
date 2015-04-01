      SUBROUTINE F08BSF(M,N,A,LDA,JPVT,TAU,WORK,RWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZGEQPF(M,N,A,LDA,JPVT,TAU,WORK,RWORK,INFO)
C
C  Purpose
C  =======
C
C  ZGEQPF computes a QR factorization with column pivoting of a
C  complex m by n matrix A: A*P = Q*R
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
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
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
C  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
C          Stores further details of
C          the orthogonal matrix Q (see A).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N)
C
C  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)
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
C  where tau is a complex scalar, and v is a complex vector with
C  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
C
C  The matrix P is represented in jpvt as follows: If
C     jpvt(j) = i
C  then the jth column of P is the ith canonical unit vector.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
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
      COMPLEX*16        A(LDA,*), TAU(*), WORK(*)
      DOUBLE PRECISION  RWORK(*)
      INTEGER           JPVT(*)
C     .. Local Scalars ..
      COMPLEX*16        AII
      DOUBLE PRECISION  TEMP, TEMP2
      INTEGER           I, ITEMP, J, MA, MN, PVT
C     .. External Functions ..
      DOUBLE PRECISION  DZNRM2
      INTEGER           IDAMAX
      EXTERNAL          DZNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASV, F08ASW, F08ASZ, F08AUZ, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCMPLX, DCONJG, MAX, MIN, SQRT
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
         CALL F06AAZ('F08BSF/ZGEQPF',-INFO)
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
               CALL ZSWAP(M,A(1,I),1,A(1,ITEMP),1)
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
         CALL F08ASZ(M,MA,A,LDA,TAU,WORK,INFO)
         IF (MA.LT.N) THEN
            CALL F08AUZ('Left','Conjugate transpose',M,N-MA,MA,A,LDA,
     *                  TAU,A(1,MA+1),LDA,WORK,INFO)
         END IF
      END IF
C
      IF (ITEMP.LT.MN) THEN
C
C        Initialize partial column norms. The first n entries of
C        work store the exact column norms.
C
         DO 40 I = ITEMP + 1, N
            RWORK(I) = DZNRM2(M-ITEMP,A(ITEMP+1,I),1)
            RWORK(N+I) = RWORK(I)
   40    CONTINUE
C
C        Compute factorization
C
         DO 80 I = ITEMP + 1, MN
C
C           Determine ith pivot column and swap if necessary
C
            PVT = (I-1) + IDAMAX(N-I+1,RWORK(I),1)
C
            IF (PVT.NE.I) THEN
               CALL ZSWAP(M,A(1,PVT),1,A(1,I),1)
               ITEMP = JPVT(PVT)
               JPVT(PVT) = JPVT(I)
               JPVT(I) = ITEMP
               RWORK(PVT) = RWORK(I)
               RWORK(N+PVT) = RWORK(N+I)
            END IF
C
C           Generate elementary reflector H(i)
C
            AII = A(I,I)
            CALL F08ASV(M-I+1,AII,A(MIN(I+1,M),I),1,TAU(I))
            A(I,I) = AII
C
            IF (I.LT.N) THEN
C
C              Apply H(i) to A(i:m,i+1:n) from the left
C
               AII = A(I,I)
               A(I,I) = DCMPLX(ONE)
               CALL F08ASW('Left',M-I+1,N-I,A(I,I),1,DCONJG(TAU(I)),
     *                     A(I,I+1),LDA,WORK)
               A(I,I) = AII
            END IF
C
C           Update partial column norms
C
            DO 60 J = I + 1, N
               IF (RWORK(J).NE.ZERO) THEN
                  TEMP = ONE - (ABS(A(I,J))/RWORK(J))**2
                  TEMP = MAX(TEMP,ZERO)
                  TEMP2 = ONE + 0.05D0*TEMP*(RWORK(J)/RWORK(N+J))**2
                  IF (TEMP2.EQ.ONE) THEN
                     IF (M-I.GT.0) THEN
                        RWORK(J) = DZNRM2(M-I,A(I+1,J),1)
                        RWORK(N+J) = RWORK(J)
                     ELSE
                        RWORK(J) = ZERO
                        RWORK(N+J) = ZERO
                     END IF
                  ELSE
                     RWORK(J) = RWORK(J)*SQRT(TEMP)
                  END IF
               END IF
   60       CONTINUE
C
   80    CONTINUE
      END IF
      RETURN
C
C     End of F08BSF (ZGEQPF)
C
      END
