      SUBROUTINE F08KSZ(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZGEBD2(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZGEBD2 reduces a complex general m by n matrix A to upper or lower
C  real bidiagonal form B by a unitary transformation: Q' * A * P = B.
C
C  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows in the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns in the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the m by n general matrix to be reduced.
C          On exit,
C          if m >= n, the diagonal and the first superdiagonal are
C            overwritten with the upper bidiagonal matrix B; the
C            elements below the diagonal, with the array TAUQ, represent
C            the unitary matrix Q as a product of elementary
C            reflectors, and the elements above the first superdiagonal,
C            with the array TAUP, represent the unitary matrix P as
C            a product of elementary reflectors;
C          if m < n, the diagonal and the first subdiagonal are
C            overwritten with the lower bidiagonal matrix B; the
C            elements below the first subdiagonal, with the array TAUQ,
C            represent the unitary matrix Q as a product of
C            elementary reflectors, and the elements above the diagonal,
C            with the array TAUP, represent the unitary matrix P as
C            a product of elementary reflectors.
C          See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
C          The diagonal elements of the bidiagonal matrix B:
C          D(i) = A(i,i).
C
C  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
C          The off-diagonal elements of the bidiagonal matrix B:
C          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
C          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
C
C  TAUQ    (output) COMPLEX*16 array dimension (min(M,N))
C          The scalar factors of the elementary reflectors which
C          represent the unitary matrix Q. See Further Details.
C
C  TAUP    (output) COMPLEX*16 array, dimension (min(M,N))
C          The scalar factors of the elementary reflectors which
C          represent the unitary matrix P. See Further Details.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (max(M,N))
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  The matrices Q and P are represented as products of elementary
C  reflectors:
C
C  If m >= n,
C
C     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
C
C  Each H(i) and G(i) has the form:
C
C     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
C
C  where tauq and taup are complex scalars, and v and u are complex
C  vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in
C  A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in
C  A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
C
C  If m < n,
C
C     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
C
C  Each H(i) and G(i) has the form:
C
C     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
C
C  where tauq and taup are complex scalars, v and u are complex vectors;
C  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
C  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
C  tauq is stored in TAUQ(i) and taup in TAUP(i).
C
C  The contents of A on exit are illustrated by the following examples:
C
C  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
C
C    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
C    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
C    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
C    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
C    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
C    (  v1  v2  v3  v4  v5 )
C
C  where d and e denote diagonal and off-diagonal elements of B, vi
C  denotes an element of the vector defining H(i), and ui an element of
C  the vector defining G(i).
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
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAUP(*), TAUQ(*), WORK(*)
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07FRY, F08ASV, F08ASW
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -4
      END IF
      IF (INFO.LT.0) THEN
         CALL F06AAZ('F08KSZ/ZGEBD2',-INFO)
         RETURN
      END IF
C
      IF (M.GE.N) THEN
C
C        Reduce to upper bidiagonal form
C
         DO 20 I = 1, N
C
C           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
C
            ALPHA = A(I,I)
            CALL F08ASV(M-I+1,ALPHA,A(MIN(I+1,M),I),1,TAUQ(I))
            D(I) = ALPHA
            A(I,I) = ONE
C
C           Apply H(i)' to A(i:m,i+1:n) from the left
C
            CALL F08ASW('Left',M-I+1,N-I,A(I,I),1,DCONJG(TAUQ(I)),
     *                  A(I,I+1),LDA,WORK)
            A(I,I) = D(I)
C
            IF (I.LT.N) THEN
C
C              Generate elementary reflector G(i) to annihilate
C              A(i,i+2:n)
C
               CALL F07FRY(N-I,A(I,I+1),LDA)
               ALPHA = A(I,I+1)
               CALL F08ASV(N-I,ALPHA,A(I,MIN(I+2,N)),LDA,TAUP(I))
               E(I) = ALPHA
               A(I,I+1) = ONE
C
C              Apply G(i) to A(i+1:m,i+1:n) from the right
C
               CALL F08ASW('Right',M-I,N-I,A(I,I+1),LDA,TAUP(I),
     *                     A(I+1,I+1),LDA,WORK)
               CALL F07FRY(N-I,A(I,I+1),LDA)
               A(I,I+1) = E(I)
            ELSE
               TAUP(I) = ZERO
            END IF
   20    CONTINUE
      ELSE
C
C        Reduce to lower bidiagonal form
C
         DO 40 I = 1, M
C
C           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
C
            CALL F07FRY(N-I+1,A(I,I),LDA)
            ALPHA = A(I,I)
            CALL F08ASV(N-I+1,ALPHA,A(I,MIN(I+1,N)),LDA,TAUP(I))
            D(I) = ALPHA
            A(I,I) = ONE
C
C           Apply G(i) to A(i+1:m,i:n) from the right
C
            CALL F08ASW('Right',M-I,N-I+1,A(I,I),LDA,TAUP(I),
     *                  A(MIN(I+1,M),I),LDA,WORK)
            CALL F07FRY(N-I+1,A(I,I),LDA)
            A(I,I) = D(I)
C
            IF (I.LT.M) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(i+2:m,i)
C
               ALPHA = A(I+1,I)
               CALL F08ASV(M-I,ALPHA,A(MIN(I+2,M),I),1,TAUQ(I))
               E(I) = ALPHA
               A(I+1,I) = ONE
C
C              Apply H(i)' to A(i+1:m,i+1:n) from the left
C
               CALL F08ASW('Left',M-I,N-I,A(I+1,I),1,DCONJG(TAUQ(I)),
     *                     A(I+1,I+1),LDA,WORK)
               A(I+1,I) = E(I)
            ELSE
               TAUQ(I) = ZERO
            END IF
   40    CONTINUE
      END IF
      RETURN
C
C     End of F08KSZ (ZGEBD2)
C
      END
