      SUBROUTINE F07ARG(M,N,A,LDA,PIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1016 (JUN 1993).
C
C  Purpose
C  =======
C
C  F07ARG computes an LU factorization of a general m-by-n matrix A
C  using partial pivoting with row interchanges.
C
C  The factorization has the form
C     A = P * L * U
C  where P is a permutation matrix, L is lower triangular (lower
C  trapezoidal if m > n), and U is upper triangular with unit diagonal
C  elements (upper trapezoidal if m < n).
C
C  This is the Level 3 BLAS version of the right-looking algorithm.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the m by n matrix to be factored.
C          On exit, the factors L and U from the factorization
C          A = P*L*U; the unit diagonal elements of U are not stored.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  PIV     (output) REAL array, dimension (M)
C          The pivot indices; for 1 <= i <= min(M,N), row i of the
C          matrix was interchanged with row PIV(i). The rest of PIV is
C          used for workspace.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, L(k,k) is exactly zero; the factorization
C               has not been completed.
C
C  This is a modified version of the LAPACK routine F07ARF/ZGETRF, in
C  which the INTEGER array IPIV has been replaced by a REAL array PIV,
C  row-equilibration is used in the choice of pivot, U has unit diagonal
C  elements, and the routine exits immediately if singularity is
C  detected.
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, RONE
      PARAMETER         (ZERO=0.0D+0,RONE=1.0D+0)
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      DOUBLE PRECISION  PIV(*)
C     .. Local Scalars ..
      INTEGER           I, IINFO, J, JB, NB
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ARH, F07ARJ, F07ZAZ, ZGEMM, ZTRSM
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DIMAG, MAX, MIN, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters.
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
         CALL F06AAZ('F07ARG       ',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Compute 2-norm of each row and store reciprocal in PIV.
C
      DO 20 I = 1, M
         PIV(I) = ZERO
   20 CONTINUE
      DO 60 J = 1, N
         DO 40 I = 1, M
            PIV(I) = PIV(I) + DBLE(A(I,J))**2 + DIMAG(A(I,J))**2
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, M
         IF (PIV(I).LE.ZERO) THEN
            INFO = I
            RETURN
         ELSE
            PIV(I) = RONE/SQRT(PIV(I))
         END IF
   80 CONTINUE
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07ARG',NB,0)
      IF (NB.LE.1) THEN
C
C        Use unblocked code.
C
         CALL F07ARH(M,N,A,LDA,PIV,INFO)
      ELSE
C
C        Use blocked code.
C
         DO 120 J = 1, MIN(M,N), NB
            JB = MIN(MIN(M,N)-J+1,NB)
C
C           Factorize diagonal and subdiagonal blocks and test for
C           exact singularity.
C
            CALL F07ARH(M-J+1,JB,A(J,J),LDA,PIV(J),IINFO)
C
            IF (IINFO.GT.0) THEN
               INFO = IINFO + J - 1
               RETURN
            END IF
C
C           Update pivot indices and apply the interchanges to columns
C           1:J-1.
C
            DO 100 I = J, MIN(M,J+JB-1)
               PIV(I) = J - 1 + PIV(I)
  100       CONTINUE
            CALL F07ARJ(J-1,A,LDA,J,J+JB-1,PIV,1)
C
            IF (J+JB.LE.N) THEN
C
C              Apply the interchanges to columns J+JB:N.
C
               CALL F07ARJ(N-J-JB+1,A(1,J+JB),LDA,J,J+JB-1,PIV,1)
C
C              Compute block row of U.
C
               CALL ZTRSM('Left','Lower','No transpose','Non-unit',JB,
     *                    N-J-JB+1,ONE,A(J,J),LDA,A(J,J+JB),LDA)
               IF (J+JB.LE.M) THEN
C
C                 Update trailing submatrix.
C
                  CALL ZGEMM('No transpose','No transpose',M-J-JB+1,
     *                       N-J-JB+1,JB,-ONE,A(J+JB,J),LDA,A(J,J+JB),
     *                       LDA,ONE,A(J+JB,J+JB),LDA)
               END IF
            END IF
  120    CONTINUE
      END IF
      RETURN
C
C     End of F07ARG
C
      END
