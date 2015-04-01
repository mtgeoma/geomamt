      SUBROUTINE F07TJZ(UPLO,DIAG,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             DTRTI2(UPLO,DIAG,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  DTRTI2 computes the inverse of a real upper or lower triangular
C  matrix A.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the matrix A is upper or lower triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  DIAG    (input) CHARACTER*1
C          Specifies whether or not the matrix A is unit triangular.
C          = 'N':  Non-unit triangular
C          = 'U':  Unit triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the triangular matrix A.  If UPLO = 'U', the
C          leading n by n upper triangular part of the array A contains
C          the upper triangular matrix, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n by n lower triangular part of the array A contains
C          the lower triangular matrix, and the strictly upper
C          triangular part of A is not referenced.  If DIAG = 'U', the
C          diagonal elements of A are also not referenced and are
C          assumed to be 1.
C          On exit, the (triangular) inverse of the original matrix, in
C          the same storage format.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the k-th diagonal element of A is zero,
C               indicating that the matrix is singular and the inverse
C               has not been computed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         DIAG, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJJ
      INTEGER           J
      LOGICAL           NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, DSCAL, DTRMV
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF ( .NOT. NOUNIT .AND. .NOT. (DIAG.EQ.'U' .OR. DIAG.EQ.'u'))
     *         THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07TJZ/DTRTI2',-INFO)
         RETURN
      END IF
C
C     Check for singularity.
C
      IF (NOUNIT) THEN
         DO 20 INFO = 1, N
            IF (A(INFO,INFO).EQ.ZERO) RETURN
   20    CONTINUE
      END IF
      INFO = 0
C
      IF (UPPER) THEN
C
C        Compute inverse of upper triangular matrix.
C
         DO 40 J = 1, N
            IF (NOUNIT) THEN
               A(J,J) = ONE/A(J,J)
               AJJ = -A(J,J)
            ELSE
               AJJ = -ONE
            END IF
C
C           Compute elements 1:j-1 of j-th column.
C
            CALL DTRMV('Upper','No transpose',DIAG,J-1,A,LDA,A(1,J),1)
            CALL DSCAL(J-1,AJJ,A(1,J),1)
   40    CONTINUE
      ELSE
C
C        Compute inverse of lower triangular matrix.
C
         DO 60 J = N, 1, -1
            IF (NOUNIT) THEN
               A(J,J) = ONE/A(J,J)
               AJJ = -A(J,J)
            ELSE
               AJJ = -ONE
            END IF
            IF (J.LT.N) THEN
C
C              Compute elements j+1:n of j-th column.
C
               CALL DTRMV('Lower','No transpose',DIAG,N-J,A(J+1,J+1),
     *                    LDA,A(J+1,J),1)
               CALL DSCAL(N-J,AJJ,A(J+1,J),1)
            END IF
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of F07TJZ (DTRTI2)
C
      END
