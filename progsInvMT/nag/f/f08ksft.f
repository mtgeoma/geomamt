      SUBROUTINE F08KSF(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZGEBRD(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZGEBRD reduces a complex general m by n matrix A to upper or lower
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
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The length of the array WORK.  LWORK >= max(1,M,N).
C          For optimum performance LWORK should be at least (M+N)*NB,
C          where NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit.
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
C  where tauq and taup are complex scalars, and v and u are complex
C  vectors; v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in
C  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in
C  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
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
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAUP(*), TAUQ(*), WORK(LWORK)
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      INTEGER           I, IINFO, IWS, J, LDWRKX, LDWRKY, MINMN, NB,
     *                  NBMIN, NX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08KSY, F08KSZ, ZGEMM
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
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
      ELSE IF (LWORK.LT.MAX(1,M,N)) THEN
         INFO = -10
      END IF
      IF (INFO.LT.0) THEN
         CALL F06AAZ('F08KSF/ZGEBRD',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      MINMN = MIN(M,N)
      IF (MINMN.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
      IWS = MAX(M,N)
      LDWRKX = M
      LDWRKY = N
C
C     Set the block size NB and the crossover point NX.
C
      CALL F07ZAZ(1,'F08KSF',NB,0)
      NB = MAX(1,NB)
C
      IF (NB.GT.1 .AND. NB.LT.MINMN) THEN
C
C        Determine when to switch from blocked to unblocked code.
C
         CALL F07ZAZ(3,'F08KSF',NX,0)
         NX = MAX(NB,NX)
         IF (NX.LT.MINMN) THEN
            IWS = (M+N)*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough work space for the optimal NB, consider using
C              a smaller block size.
C
               CALL F07ZAZ(2,'F08KSF',NBMIN,0)
               IF (LWORK.GE.(M+N)*NBMIN) THEN
                  NB = LWORK/(M+N)
               ELSE
                  NB = 1
                  NX = MINMN
               END IF
            END IF
         END IF
      ELSE
         NX = MINMN
      END IF
C
      DO 60 I = 1, MINMN - NX, NB
C
C        Reduce rows and columns i:i+ib-1 to bidiagonal form and return
C        the matrices X and Y which are needed to update the unreduced
C        part of the matrix
C
         CALL F08KSY(M-I+1,N-I+1,NB,A(I,I),LDA,D(I),E(I),TAUQ(I),TAUP(I)
     *               ,WORK,LDWRKX,WORK(LDWRKX*NB+1),LDWRKY)
C
C        Update the trailing submatrix A(i+ib:m,i+ib:n), using
C        an update of the form  A := A - V*Y' - X*U'
C
         CALL ZGEMM('No transpose','Conjugate transpose',M-I-NB+1,
     *              N-I-NB+1,NB,-ONE,A(I+NB,I),LDA,WORK(LDWRKX*NB+NB+1),
     *              LDWRKY,ONE,A(I+NB,I+NB),LDA)
         CALL ZGEMM('No transpose','No transpose',M-I-NB+1,N-I-NB+1,NB,
     *              -ONE,WORK(NB+1),LDWRKX,A(I,I+NB),LDA,ONE,
     *              A(I+NB,I+NB),LDA)
C
C        Copy diagonal and off-diagonal elements of B back into A
C
         IF (M.GE.N) THEN
            DO 20 J = I, I + NB - 1
               A(J,J) = D(J)
               A(J,J+1) = E(J)
   20       CONTINUE
         ELSE
            DO 40 J = I, I + NB - 1
               A(J,J) = D(J)
               A(J+1,J) = E(J)
   40       CONTINUE
         END IF
   60 CONTINUE
C
C     Use unblocked code to reduce the remainder of the matrix
C
      CALL F08KSZ(M-I+1,N-I+1,A(I,I),LDA,D(I),E(I),TAUQ(I),TAUP(I),WORK,
     *            IINFO)
      WORK(1) = IWS
      RETURN
C
C     End of F08KSF (ZGEBRD)
C
      END
