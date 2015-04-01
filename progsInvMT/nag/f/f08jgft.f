      SUBROUTINE F08JGF(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DPTEQR(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
C
C  Purpose
C  =======
C
C  DPTEQR computes all eigenvalues and, optionally, eigenvectors of a
C  symmetric positive definite tridiagonal matrix by first factoring the
C  matrix using F08JGZ and then calling F08MEF to compute the singular
C  values of the bidiagonal factor.  This routine computes the
C  eigenvalues of the positive definite tridiagonal matrix to high
C  relative accuracy.  This means that if the eigenvalues range over
C  many orders of magnitude in size, then the small eigenvalues and
C  corresponding eigenvectors will be computed more accurately than, for
C  example, with the standard QR method.  The eigenvectors of a full or
C  band symmetric matrix can also be found if DSYTRD or DSPTRD or DSBTRD
C  has been used to reduce this matrix to tridiagonal form.
C  (The reduction to tridiagonal form, however, may preclude the
C  possibility of obtaining high relative accuracy in the small
C  eigenvalues of the original matrix, if these eigenvalues range over
C  many orders of magnitude.)
C
C  Arguments
C  =========
C
C  COMPZ   (input) CHARACTER*1
C          Specifies whether eigenvectors are to be computed
C          as follows
C
C             COMPZ = 'N' or 'n'   Compute eigenvalues only.
C
C             COMPZ = 'V' or 'v'   Compute eigenvectors of original
C                                  symmetric matrix also.
C                                  Array Z contains the orthogonal
C                                  matrix used to reduce the original
C                                  matrix to tridiagonal form.
C
C             COMPZ = 'I' or 'i'   Compute eigenvectors of
C                                  tridiagonal matrix also.
C
C  N       (input) INTEGER
C          The number of rows and columns in the matrix.  N >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, D contains the diagonal elements of the
C          tridiagonal matrix.
C          On normal exit, D contains the eigenvalues, in descending
C          order.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, E contains the subdiagonal elements of the
C          tridiagonal matrix in positions 1 through N-1.
C          E(N) is arbitrary.
C          On exit, E has been destroyed.
C
C  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
C          If  COMPZ = 'V' or 'v', then:
C          On entry, Z contains the orthogonal matrix used in the
C          reduction to tridiagonal form.
C          If  COMPZ = 'V' or 'v' or 'I' or 'i', then:
C          On exit, Z contains the orthonormal eigenvectors of the
C          symmetric tridiagonal (or full) matrix.  If an error exit
C          is made, Z contains the eigenvectors associated with the
C          stored eigenvalues.
C
C          If  COMPZ = 'N' or 'n', then Z is not referenced.
C
C  LDZ     (input) INTEGER
C          The leading dimension of the array Z.  If eigenvectors are
C          desired, then  LDZ >= max( 1, N ).  In any case, LDZ >= 1.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,4*N-4))
C          Workspace used in computing eigenvectors.
C          If  COMPZ = 'N' or 'n', then WORK is not referenced.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit.
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C          > 0:  if INFO = +i, 1 <= i <= N, the Cholesky factorization
C                              of the matrix could not be performed
C                              because the i-th principal minor was not
C                              positive definite.
C                if INFO = N+i, 1 <= i <= N, the i-th singular value
C                              of the bidiagonal factor failed to
C                              converge.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDZ, N
      CHARACTER         COMPZ
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E(*), WORK(*), Z(LDZ,*)
C     .. Local Scalars ..
      INTEGER           I, ICOMPZ, II, J, NRU
C     .. Local Arrays ..
      DOUBLE PRECISION  C(1,1), VT(1,1)
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DBDSQR, DSCAL, F06AAZ, F06QHF, F08JGZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
C
      IF ((COMPZ.EQ.'N' .OR. COMPZ.EQ.'n')) THEN
         ICOMPZ = 0
      ELSE IF ((COMPZ.EQ.'V' .OR. COMPZ.EQ.'v')) THEN
         ICOMPZ = 1
      ELSE IF ((COMPZ.EQ.'I' .OR. COMPZ.EQ.'i')) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF (ICOMPZ.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF ((LDZ.LT.1) .OR. (ICOMPZ.GT.0 .AND. LDZ.LT.MAX(1,N))) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08JGF/DPTEQR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF (N.EQ.1) THEN
         IF (ICOMPZ.GT.0) Z(1,1) = ONE
         RETURN
      END IF
      IF (ICOMPZ.EQ.2) CALL F06QHF('General',N,N,ZERO,ONE,Z,LDZ)
C
C     Call F08JGZ to factor the matrix.
C
      CALL F08JGZ(N,D,E,INFO)
      IF (INFO.NE.0) RETURN
      DO 20 I = 1, N
         D(I) = SQRT(D(I))
   20 CONTINUE
      DO 40 I = 1, N - 1
         E(I) = E(I)*D(I)
   40 CONTINUE
C
C     Call F08MEF to compute the singular values/vectors of the
C     bidiagonal factor.
C
      IF (ICOMPZ.GT.0) THEN
         NRU = N
      ELSE
         NRU = 0
      END IF
      CALL DBDSQR('Lower',N,0,NRU,0,D,E,VT,1,Z,LDZ,C,1,WORK,INFO)
C
C     Square the singular values.
C
      IF (INFO.EQ.0) THEN
         DO 60 I = 1, N
            D(I) = D(I)*D(I)
   60    CONTINUE
      ELSE
         INFO = N + INFO
      END IF
C
      IF (ICOMPZ.GT.0 .AND. INFO.EQ.0) THEN
C
C        Normalize eigenvectors so that element of largest absolute
C        value is positive
C
         DO 80 J = 1, N
            II = IDAMAX(N,Z(1,J),1)
            IF (Z(II,J).LT.ZERO) CALL DSCAL(N,-ONE,Z(1,J),1)
   80    CONTINUE
      END IF
C
      RETURN
C
C     End of F08JGF (DPTEQR)
C
      END
