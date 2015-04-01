      SUBROUTINE F08DVZ(M,P,N,A,LDA,TAUA,B,LDB,TAUB,WORK,LWORK,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  ZGGRQF computes a generalized RQ factorization of an M-by-N matrix A
C  and a P-by-N matrix B:
C
C              A = R*Q,        B = Z*T*Q,
C
C  where Q is an N-by-N unitary matrix, Z is a P-by-P unitary
C  matrix, and R and T assume one of the forms:
C
C  if M <= N,  R = ( 0  R12 ) M,   or if M > N,  R = ( R11 ) M-N,
C                   N-M  M                           ( R21 ) N
C                                                       N
C
C  where R12 or R21 is upper triangular, and
C
C  if P >= N,  T = ( T11 ) N  ,   or if P < N,  T = ( T11  T12 ) P,
C                  (  0  ) P-N                         P   N-P
C                     N
C
C  where T11 is upper triangular.
C
C  In particular, if B is square and nonsingular, the GRQ factorization
C  of A and B implicitly gives the RQ factorization of A*inv(B):
C
C               A*inv(B) = (R*inv(T))*Z'
C
C  where inv(B) denotes the inverse of the matrix B, and Z' denotes the
C  conjugate transpose of the matrix Z.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  P       (input) INTEGER
C          The number of rows of the matrix B.  P >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrices A and B. N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the M-by-N matrix A.
C          On exit, if M <= N, the upper triangle of the subarray
C          A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R;
C          if M > N, the elements on and above the (M-N)-th subdiagonal
C          contain the M-by-N upper trapezoidal matrix R; the remaining
C          elements, with the array TAUA, represent the unitary
C          matrix Q as a product of elementary reflectors (see Further
C          Details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,M).
C
C  TAUA    (output) COMPLEX*16 array, dimension (min(M,N))
C          The scalar factors of the elementary reflectors which
C          represent the unitary matrix Q (see Further Details).
C
C  B       (input/output) COMPLEX*16 array, dimension (LDB,N)
C          On entry, the P-by-N matrix B.
C          On exit, the elements on and above the diagonal of the array
C          contain the min(P,N)-by-N upper trapezoidal matrix T (T is
C          upper triangular if P >= N); the elements below the diagonal,
C          with the array TAUB, represent the unitary matrix Z as a
C          product of elementary reflectors (see Further Details).
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B. LDB >= max(1,P).
C
C  TAUB    (output) COMPLEX*16 array, dimension (min(P,N))
C          The scalar factors of the elementary reflectors which
C          represent the unitary matrix Z (see Further Details).
C
C  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,N,M,P).
C          For optimum performance LWORK >= max(M*NB1,N*NB2,P*NB3),
C          where NB1 is the optimal blocksize for the RQ factorization
C          of an M-by-N matrix, NB2 is the optimal blocksize for the
C          QR factorization of a P-by-N matrix, and NB3 is the optimal
C          blocksize for a call of ZUNMRQ.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO=-i, the i-th argument had an illegal value.
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
C     H(i) = I - taua * v * v'
C
C  where taua is a complex scalar, and v is a complex vector with
C  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
C  A(m-k+i,1:n-k+i-1), and taua in TAUA(i).
C  To form Q explicitly, use LAPACK subroutine ZUNGRQ.
C  To use Q to update another matrix, use LAPACK subroutine ZUNMRQ.
C
C  The matrix Z is represented as a product of elementary reflectors
C
C     Z = H(1) H(2) . . . H(k), where k = min(p,n).
C
C  Each H(i) has the form
C
C     H(i) = I - taub * v * v'
C
C  where taub is a complex scalar, and v is a complex vector with
C  v(1:i-1) = 0 and v(i) = 1; v(i+1:p) is stored on exit in B(i+1:p,i),
C  and taub in TAUB(i).
C  To form Z explicitly, use LAPACK subroutine ZUNGQR.
C  To use Z to update another matrix, use LAPACK subroutine ZUNMQR.
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LWORK, M, N, P
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), TAUA(*), TAUB(*),
     *                  WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           LOPT
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08CVZ, F08CXZ, ZGEQRF
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (P.LT.0) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      ELSE IF (LDB.LT.MAX(1,P)) THEN
         INFO = -8
      ELSE IF (LWORK.LT.MAX(1,M,P,N)) THEN
         INFO = -11
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08DVZ',-INFO)
         RETURN
      END IF
C
C     RQ factorization of M-by-N matrix A: A = R*Q
C
      CALL F08CVZ(M,N,A,LDA,TAUA,WORK,LWORK,INFO)
      LOPT = WORK(1)
C
C     Update B := B*Q'
C
      CALL F08CXZ('Right','Conjugate Transpose',P,N,MIN(M,N),
     *            A(MAX(1,M-N+1),1),LDA,TAUA,B,LDB,WORK,LWORK,INFO)
      LOPT = MAX(LOPT,INT(WORK(1)))
C
C     QR factorization of P-by-N matrix B: B = Z*T
C
      CALL ZGEQRF(P,N,B,LDB,TAUB,WORK,LWORK,INFO)
      WORK(1) = MAX(LOPT,INT(WORK(1)))
C
      RETURN
C
C     End of F08DVZ (ZGGRQF)
C
      END
