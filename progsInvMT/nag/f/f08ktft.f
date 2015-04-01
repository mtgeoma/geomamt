      SUBROUTINE F08KTF(VECT,M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZUNGBR(VECT,M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNGBR generates one of the matrices Q or P' determined by ZGEBRD
C  when reducing a complex matrix A to bidiagonal form: A = Q * B * P'.
C  Q and P' are defined as products of elementary reflectors H(i) or
C  G(i) respectively..
C
C  If VECT = 'Q', A is assumed to have been an m-by-k matrix, and Q
C  is of order m:
C  if m >= k, Q = H(1) H(2) . . . H(k) and ZUNGBR returns the first n
C  columns of Q, where m >= n >= k;
C  if m < k, Q = H(1) H(2) . . . H(m-1) and ZUNGBR returns Q as an
C  m-by-m matrix.
C
C  If VECT = 'P', A is assumed to have been a k-by-n matrix, and P'
C  is of order n:
C  if k < n, P' = G(k) . . . G(2) G(1) and ZUNGBR returns the first m
C  rows of P', where n >= m >= k;
C  if k >= n, P' = G(n-1) . . . G(2) G(1) and ZUNGBR returns P' as an
C  n-by-n matrix.
C
C  Arguments
C  =========
C
C  VECT    (input) CHARACTER*1
C          Specifies whether the matrix Q or the matrix P' is required,
C          as defined in the transformation applied by ZGEBRD:
C          = 'Q': generate Q;
C          = 'P': generate P'.
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q or P' to be returned.
C          M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q or P' to be returned.
C          N >= 0.
C          If VECT = 'Q', M >= N >= min(M,K);
C          if VECT = 'P', N >= M >= min(N,K).
C
C  K       (input) INTEGER
C          If VECT = 'Q', the number of columns in the original m-by-k
C          matrix reduced by ZGEBRD.
C          If VECT = 'P', the number of rows in the original k-by-n
C          matrix reduced by ZGEBRD.
C          K >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the vectors which define the elementary reflectors,
C          as returned by ZGEBRD.
C          On exit, the m-by-n matrix Q or P'.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= M.
C
C  TAU     (input) COMPLEX*16 array, dimension
C                                (min(M,K)) if VECT = 'Q'
C                                (min(N,K)) if VECT = 'P'
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i) or G(i), which determines Q or P', as returned
C          by ZGEBRD in its array argument TAUQ or TAUP.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
C          For optimum performance LWORK should be at least min(M,N)*NB,
C          where NB is the optimal blocksize.
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
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, LWORK, M, N
      CHARACTER         VECT
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IINFO, J
      LOGICAL           WANTQ
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZUNGLQ, ZUNGQR
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      WANTQ = (VECT.EQ.'Q' .OR. VECT.EQ.'q')
      IF ( .NOT. WANTQ .AND. .NOT. (VECT.EQ.'P' .OR. VECT.EQ.'p')) THEN
         INFO = -1
      ELSE IF (M.LT.0) THEN
         INFO = -2
      ELSE IF (N.LT.0 .OR. (WANTQ .AND. (N.GT.M .OR. N.LT.MIN(M,K)))
     *         .OR. ( .NOT. WANTQ .AND. (M.GT.N .OR. M.LT.MIN(N,K))))
     *         THEN
         INFO = -3
      ELSE IF (K.LT.0) THEN
         INFO = -4
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -6
      ELSE IF (LWORK.LT.MAX(1,MIN(M,N))) THEN
         INFO = -9
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08KTF/ZUNGBR',-INFO)
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
      IF (WANTQ) THEN
C
C        Form Q, determined by a call to ZGEBRD to reduce an m-by-k
C        matrix
C
         IF (M.GE.K) THEN
C
C           If m >= k, assume m >= n >= k
C
            CALL ZUNGQR(M,N,K,A,LDA,TAU,WORK,LWORK,IINFO)
C
         ELSE
C
C           If m < k, assume m = n
C
C           Shift the vectors which define the elementary reflectors one
C           column to the right, and set the first row and column of Q
C           to those of the unit matrix
C
            DO 40 J = M, 2, -1
               A(1,J) = ZERO
               DO 20 I = J + 1, M
                  A(I,J) = A(I,J-1)
   20          CONTINUE
   40       CONTINUE
            A(1,1) = ONE
            DO 60 I = 2, M
               A(I,1) = ZERO
   60       CONTINUE
            IF (M.GT.1) THEN
C
C              Form Q(2:m,2:m)
C
               CALL ZUNGQR(M-1,M-1,M-1,A(2,2),LDA,TAU,WORK,LWORK,IINFO)
            ELSE
               WORK(1) = 1
            END IF
         END IF
      ELSE
C
C        Form P', determined by a call to ZGEBRD to reduce a k-by-n
C        matrix
C
         IF (K.LT.N) THEN
C
C           If k < n, assume k <= m <= n
C
            CALL ZUNGLQ(M,N,K,A,LDA,TAU,WORK,LWORK,IINFO)
C
         ELSE
C
C           If k >= n, assume m = n
C
C           Shift the vectors which define the elementary reflectors one
C           row downward, and set the first row and column of P' to
C           those of the unit matrix
C
            A(1,1) = ONE
            DO 80 I = 2, N
               A(I,1) = ZERO
   80       CONTINUE
            DO 120 J = 2, N
               DO 100 I = J - 1, 2, -1
                  A(I,J) = A(I-1,J)
  100          CONTINUE
               A(1,J) = ZERO
  120       CONTINUE
            IF (N.GT.1) THEN
C
C              Form P'(2:n,2:n)
C
               CALL ZUNGLQ(N-1,N-1,N-1,A(2,2),LDA,TAU,WORK,LWORK,IINFO)
            ELSE
               WORK(1) = 1
            END IF
         END IF
      END IF
      RETURN
C
C     End of F08KTF (ZUNGBR)
C
      END
