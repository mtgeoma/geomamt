      SUBROUTINE F08NVF(JOB,N,A,LDA,ILO,IHI,SCALE,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZGEBAL(JOB,N,A,LDA,ILO,IHI,SCALE,INFO)
C
C  Purpose
C  =======
C
C  ZGEBAL does permutations to isolate eigenvalues and/or balances a
C  complex general matrix.  Balancing may reduce the 1-norm of the
C  matrix and can improve the accuracy of the computed eigenvalues
C  and/or eigenvectors.
C
C  Arguments
C  =========
C
C  JOB     (input) CHARACTER*1
C          Specifies the operations to be done:
C          = 'N', do nothing with the input matrix, but set
C                 ILO = 1, IHI = N, SCALE(I) = ONE.
C          = 'P', only permute to isolate eigenvalues.
C          = 'S', only balance the input matrix.
C          = 'B', both permute and balance.
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, A contains the input matrix to be balanced.
C          On exit,  A contains the permuted and/or balanced matrix.
C          See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the matrix A. LDA >= max(1,N).
C
C  ILO     (output) INTEGER
C  IHI     (output) INTEGER
C          On exit, ILO and IHI are two integers such that A(I,J)
C          is equal to zero if
C             (1) I is greater than J and
C             (2) J = 1,...,ILO-1 or I = IHI+1,...,N.
C
C  SCALE   (output) DOUBLE PRECISION array, dimension (N)
C          On exit, SCALE contains information about the permutations
C          and scaling factors used.  If P(J) is the index of the
C          row and column interchanged with row and column J and D(J,J)
C          is the scaling factor used to balance row and column J of the
C          submatrix of A, then
C          SCALE(J) = P(J),    for J = 1,...,ILO-1
C                   = D(J,J),      J = ILO,...,IHI
C                   = P(J)         J = IHI+1,...,N.
C          The order in which the interchanges are made is N to IHI+1,
C          then 1 to ILO-1.
C
C  INFO    (output) INTEGER
C          = 0:  normal return.
C          < 0:  if INFO = -k, the k-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  The permutations consist of row and column interchanges which put
C  the matrix in the form
C
C             ( T1   X   Y  )
C     P A P = (  0   B   Z  )
C             (  0   0   T2 )
C
C  where T1 and T2 are upper triangular matrices whose eigenvalues lie
C  along the diagonal.  The column indices ILO and IHI mark the starting
C  and ending columns of the submatrix B. Balancing consists of applying
C  a diagonal similarity transformation inv(D) * B * D to make the
C  1-norms of each row of B and its corresponding column nearly equal.
C  The output matrix is
C
C     ( T1     X*D          Y    )
C     (  0  inv(D)*B*D  inv(D)*Z ).
C     (  0      0           T2   )
C
C  Information about the permutations P and the diagonal matrix D is
C  returned in the vector SCALE.
C
C  This subroutine is based on the EISPACK routine CBAL.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  SCLFAC, B2
      PARAMETER         (SCLFAC=1.0D+1,B2=SCLFAC*SCLFAC)
      DOUBLE PRECISION  FACTOR
      PARAMETER         (FACTOR=0.95D+0)
C     .. Scalar Arguments ..
      INTEGER           IHI, ILO, INFO, LDA, N
      CHARACTER         JOB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      DOUBLE PRECISION  SCALE(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM
      DOUBLE PRECISION  C, F, G, R, S
      INTEGER           I, IEXC, J, K, L, M
      LOGICAL           NOCONV
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZDSCAL, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DIMAG, MAX
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(CDUM) = ABS(DBLE(CDUM)) + ABS(DIMAG(CDUM))
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF ( .NOT. (JOB.EQ.'N' .OR. JOB.EQ.'n')
     *    .AND. .NOT. (JOB.EQ.'P' .OR. JOB.EQ.'p')
     *    .AND. .NOT. (JOB.EQ.'S' .OR. JOB.EQ.'s')
     *    .AND. .NOT. (JOB.EQ.'B' .OR. JOB.EQ.'b')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08NVF/ZGEBAL',-INFO)
         RETURN
      END IF
C
      K = 1
      L = N
C
      IF (N.EQ.0) GO TO 420
C
      IF ((JOB.EQ.'N' .OR. JOB.EQ.'n')) THEN
         DO 20 I = 1, N
            SCALE(I) = ONE
   20    CONTINUE
         GO TO 420
      END IF
C
      IF ((JOB.EQ.'S' .OR. JOB.EQ.'s')) GO TO 240
C
C     Permutation to isolate eigenvalues if possible
C
      GO TO 100
C
C     Row and column exchange.
C
   40 CONTINUE
      SCALE(M) = J
      IF (J.EQ.M) GO TO 60
C
      CALL ZSWAP(L,A(1,J),1,A(1,M),1)
      CALL ZSWAP(N-K+1,A(J,K),LDA,A(M,K),LDA)
C
   60 CONTINUE
      GO TO (80,160) IEXC
C
C     Search for rows isolating an eigenvalue and push them down.
C
   80 CONTINUE
      IF (L.EQ.1) GO TO 420
      L = L - 1
C
  100 CONTINUE
      DO 140 J = L, 1, -1
C
         DO 120 I = 1, L
            IF (I.EQ.J) GO TO 120
            IF (DBLE(A(J,I)).NE.ZERO .OR. DIMAG(A(J,I)).NE.ZERO)
     *          GO TO 140
  120    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 40
  140 CONTINUE
C
      GO TO 180
C
C     Search for columns isolating an eigenvalue and push them left.
C
  160 CONTINUE
      K = K + 1
C
  180 CONTINUE
      DO 220 J = K, L
C
         DO 200 I = K, L
            IF (I.EQ.J) GO TO 200
            IF (DBLE(A(I,J)).NE.ZERO .OR. DIMAG(A(I,J)).NE.ZERO)
     *          GO TO 220
  200    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 40
  220 CONTINUE
C
  240 CONTINUE
      DO 260 I = K, L
         SCALE(I) = ONE
  260 CONTINUE
C
      IF ((JOB.EQ.'P' .OR. JOB.EQ.'p')) GO TO 420
C
C     Balance the submatrix in rows K to L.
C
C     Iterative loop for norm reduction
C
  280 CONTINUE
      NOCONV = .FALSE.
C
      DO 400 I = K, L
         C = ZERO
         R = ZERO
C
         DO 300 J = K, L
            IF (J.EQ.I) GO TO 300
            C = C + CABS1(A(J,I))
            R = R + CABS1(A(I,J))
  300    CONTINUE
C
C        Guard against zero C or R due to underflow.
C
         IF (C.EQ.ZERO .OR. R.EQ.ZERO) GO TO 400
         G = R/SCLFAC
         F = ONE
         S = C + R
  320    CONTINUE
         IF (C.GE.G) GO TO 340
         F = F*SCLFAC
         C = C*B2
         GO TO 320
C
  340    CONTINUE
         G = R*SCLFAC
  360    CONTINUE
         IF (C.LT.G) GO TO 380
         F = F/SCLFAC
         C = C/B2
         GO TO 360
C
C        Now balance.
C
  380    CONTINUE
         IF ((C+R)/F.GE.FACTOR*S) GO TO 400
         G = ONE/F
         SCALE(I) = SCALE(I)*F
         NOCONV = .TRUE.
C
         CALL ZDSCAL(N-K+1,G,A(I,K),LDA)
         CALL ZDSCAL(L,F,A(1,I),1)
C
  400 CONTINUE
C
      IF (NOCONV) GO TO 280
C
  420 CONTINUE
      ILO = K
      IHI = L
C
      RETURN
C
C     End of F08NVF (ZGEBAL)
C
      END
