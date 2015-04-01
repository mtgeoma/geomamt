      SUBROUTINE F08DSZ(N,M,P,A,LDA,TAUA,B,LDB,TAUB,WORK,LWORK,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  ZGGQRF computes a generalized QR factorization of an N-by-M matrix A
C  and an N-by-P matrix B:
C
C              A = Q*R,        B = Q*T*Z,
C
C  where Q is an N-by-N unitary matrix, Z is a P-by-P unitary matrix,
C  and R and T assume one of the forms:
C
C  if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N,
C                  (  0  ) N-M                         N   M-N
C                     M
C
C  where R11 is upper triangular, and
C
C  if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P,
C                   P-N  N                           ( T21 ) P
C                                                       P
C
C  where T12 or T21 is upper triangular.
C
C  In particular, if B is square and nonsingular, the GQR factorization
C  of A and B implicitly gives the QR factorization of inv(B)*A:
C
C               inv(B)*A = Z'*(inv(T)*R)
C
C  where inv(B) denotes the inverse of the matrix B, and Z' denotes the
C  conjugate transpose of matrix Z.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of rows of the matrices A and B. N >= 0.
C
C  M       (input) INTEGER
C          The number of columns of the matrix A.  M >= 0.
C
C  P       (input) INTEGER
C          The number of columns of the matrix B.  P >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,M)
C          On entry, the N-by-M matrix A.
C          On exit, the elements on and above the diagonal of the array
C          contain the min(N,M)-by-M upper trapezoidal matrix R (R is
C          upper triangular if N >= M); the elements below the diagonal,
C          with the array TAUA, represent the unitary matrix Q as a
C          product of min(N,M) elementary reflectors (see Further
C          Details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,N).
C
C  TAUA    (output) COMPLEX*16 array, dimension (min(N,M))
C          The scalar factors of the elementary reflectors which
C          represent the unitary matrix Q (see Further Details).
C
C  B       (input/output) COMPLEX*16 array, dimension (LDB,P)
C          On entry, the N-by-P matrix B.
C          On exit, if N <= P, the upper triangle of the subarray
C          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;
C          if N > P, the elements on and above the (N-P)-th subdiagonal
C          contain the N-by-P upper trapezoidal matrix T; the remaining
C          elements, with the array TAUB, represent the unitary
C          matrix Z as a product of elementary reflectors (see Further
C          Details).
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B. LDB >= max(1,N).
C
C  TAUB    (output) COMPLEX*16 array, dimension (min(N,P))
C          The scalar factors of the elementary reflectors which
C          represent the unitary matrix Z (see Further Details).
C
C  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,N,M,P).
C          For optimum performance LWORK >= max(M*NB1,N*NB2,P*NB3),
C          where NB1 is the optimal blocksize for the QR factorization
C          of an N-by-M matrix, NB2 is the optimal blocksize for the
C          RQ factorization of an N-by-P matrix, and NB3 is the optimal
C          blocksize for a call of ZUNMQR.
C
C  INFO    (output) INTEGER
C           = 0:  successful exit
C           < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  The matrix Q is represented as a product of elementary reflectors
C
C     Q = H(1) H(2) . . . H(k), where k = min(n,m).
C
C  Each H(i) has the form
C
C     H(i) = I - taua * v * v'
C
C  where taua is a complex scalar, and v is a complex vector with
C  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
C  and taua in TAUA(i).
C  To form Q explicitly, use LAPACK subroutine ZUNGQR.
C  To use Q to update another matrix, use LAPACK subroutine ZUNMQR.
C
C  The matrix Z is represented as a product of elementary reflectors
C
C     Z = H(1) H(2) . . . H(k), where k = min(n,p).
C
C  Each H(i) has the form
C
C     H(i) = I - taub * v * v'
C
C  where taub is a complex scalar, and v is a complex vector with
C  v(p-k+i+1:p) = 0 and v(p-k+i) = 1; v(1:p-k+i-1) is stored on exit in
C  B(n-k+i,1:p-k+i-1), and taub in TAUB(i).
C  To form Z explicitly, use LAPACK subroutine ZUNGRQ.
C  To use Z to update another matrix, use LAPACK subroutine ZUNMRQ.
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
      EXTERNAL          F06AAZ, F08CVZ, ZGEQRF, ZUNMQR
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF (N.LT.0) THEN
         INFO = -1
      ELSE IF (M.LT.0) THEN
         INFO = -2
      ELSE IF (P.LT.0) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -8
      ELSE IF (LWORK.LT.MAX(1,N,M,P)) THEN
         INFO = -11
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08DSZ',-INFO)
         RETURN
      END IF
C
C     QR factorization of N-by-M matrix A: A = Q*R
C
      CALL ZGEQRF(N,M,A,LDA,TAUA,WORK,LWORK,INFO)
      LOPT = WORK(1)
C
C     Update B := Q'*B.
C
      CALL ZUNMQR('Left','Conjugate Transpose',N,P,MIN(N,M),A,LDA,TAUA,
     *            B,LDB,WORK,LWORK,INFO)
      LOPT = MAX(LOPT,INT(WORK(1)))
C
C     RQ factorization of N-by-P matrix B: B = T*Z.
C
      CALL F08CVZ(N,P,B,LDB,TAUB,WORK,LWORK,INFO)
      WORK(1) = MAX(LOPT,INT(WORK(1)))
C
      RETURN
C
C     End of F08DSZ (ZGGQRF)
C
      END
