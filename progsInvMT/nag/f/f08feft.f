      SUBROUTINE F08FEF(UPLO,N,A,LDA,D,E,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1648 (JUN 1995).
C     .. Entry Points ..
      ENTRY             DSYTRD(UPLO,N,A,LDA,D,E,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  DSYTRD reduces a real symmetric matrix A to symmetric tridiagonal
C  form T by an orthogonal similarity transformation: Q' * A * Q = T.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit, if UPLO = 'U', the diagonal and first superdiagonal
C          of A are overwritten by the corresponding elements of the
C          tridiagonal matrix T, and the elements above the first
C          superdiagonal, with the array TAU, represent the orthogonal
C          matrix Q as a product of elementary reflectors; if UPLO
C          = 'L', the diagonal and first subdiagonal of A are over-
C          written by the corresponding elements of the tridiagonal
C          matrix T, and the elements below the first subdiagonal, with
C          the array TAU, represent the orthogonal matrix Q as a product
C          of elementary reflectors. See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  D       (output) DOUBLE PRECISION array, dimension (N)
C          The diagonal elements of the tridiagonal matrix T:
C          D(i) = A(i,i).
C
C  E       (output) DOUBLE PRECISION array, dimension (N-1)
C          The off-diagonal elements of the tridiagonal matrix T:
C          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
C
C  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.  LWORK >= 1.
C          For optimum performance LWORK should be at least N*NB,
C          where NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit.
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(n-1) . . . H(2) H(1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with
C  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
C  A(1:i-1,i+1), and tau in TAU(i).
C
C  If UPLO = 'L', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(1) H(2) . . . H(n-1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with
C  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
C  and tau in TAU(i).
C
C  The contents of A on exit are illustrated by the following examples
C  with n = 5:
C
C  if UPLO = 'U':                       if UPLO = 'L':
C
C    (  d   e   v2  v3  v4 )              (  d                  )
C    (      d   e   v3  v4 )              (  e   d              )
C    (          d   e   v4 )              (  v1  e   d          )
C    (              d   e  )              (  v1  v2  e   d      )
C    (                  d  )              (  v1  v2  v3  e   d  )
C
C  where d and e denote diagonal and off-diagonal elements of T, and vi
C  denotes an element of the vector defining H(i).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), D(*), E(*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IINFO, IWS, J, KK, LDWORK, NB, NBMIN, NX
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          DSYR2K, F06AAZ, F07ZAZ, F08FEY, F08FEZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.1) THEN
         INFO = -9
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08FEF/DSYTRD',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size.
C
      CALL F07ZAZ(1,'F08FEF',NB,0)
      NX = N
      IWS = 1
      IF (NB.GT.1 .AND. NB.LT.N) THEN
C
C        Determine when to cross over from blocked to unblocked code
C        (last block is always handled by unblocked code).
C
         CALL F07ZAZ(3,'F08FEF',NX,0)
         NX = MAX(NB,NX)
         IF (NX.LT.N) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough workspace to use optimal NB:  determine the
C              minimum value of NB, and reduce NB or force use of
C              unblocked code by setting NX = N.
C
               NB = MAX(LWORK/LDWORK,1)
               CALL F07ZAZ(2,'F08FEF',NBMIN,0)
               IF (NB.LT.NBMIN) NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
C
      IF (UPPER) THEN
C
C        Reduce the upper triangle of A.
C        Columns 1:kk are handled by the unblocked method.
C
         KK = N - ((N-NX+NB-1)/NB)*NB
         DO 40 I = N - NB + 1, KK + 1, -NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL F08FEY(UPLO,I+NB-1,NB,A,LDA,E,TAU,WORK,LDWORK)
C
C           Update the unreduced submatrix A(1:i-1,1:i-1), using an
C           update of the form:  A := A - V*W' - W*V'
C
            CALL DSYR2K(UPLO,'No transpose',I-1,NB,-ONE,A(1,I),LDA,WORK,
     *                  LDWORK,ONE,A,LDA)
C
C           Copy superdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 20 J = I, I + NB - 1
               A(J-1,J) = E(J-1)
               D(J) = A(J,J)
   20       CONTINUE
   40    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL F08FEZ(UPLO,KK,A,LDA,D,E,TAU,IINFO)
      ELSE
C
C        Reduce the lower triangle of A
C
         DO 80 I = 1, N - NX, NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL F08FEY(UPLO,N-I+1,NB,A(I,I),LDA,E(I),TAU(I),WORK,
     *                  LDWORK)
C
C           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
C           an update of the form:  A := A - V*W' - W*V'
C
            CALL DSYR2K(UPLO,'No transpose',N-I-NB+1,NB,-ONE,A(I+NB,I),
     *                  LDA,WORK(NB+1),LDWORK,ONE,A(I+NB,I+NB),LDA)
C
C           Copy subdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 60 J = I, I + NB - 1
               A(J+1,J) = E(J)
               D(J) = A(J,J)
   60       CONTINUE
   80    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL F08FEZ(UPLO,N-I+1,A(I,I),LDA,D(I),E(I),TAU(I),IINFO)
      END IF
C
      WORK(1) = IWS
      RETURN
C
C     End of F08FEF (DSYTRD)
C
      END
