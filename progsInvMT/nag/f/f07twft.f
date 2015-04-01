      SUBROUTINE F07TWF(UPLO,DIAG,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  ZTRTRI computes the inverse of a complex upper or lower triangular
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
C  A       (input/output) COMPLEX array, dimension (LDA,N)
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
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         DIAG, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
C     .. Local Scalars ..
      INTEGER           J, JB, NB, NN
      LOGICAL           NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          ZTRMM, ZTRSM, F06AAZ, F07TWZ, F07ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
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
         CALL F06AAZ('F07TWF/ZTRTRI',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
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
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07TWF',NB,0)
      IF (NB.LE.1) NB = MAX(1,N)
C
      IF (UPPER) THEN
C
C        Compute inverse of upper triangular matrix
C
         DO 40 J = 1, N, NB
            JB = MIN(NB,N-J+1)
C
C           Compute rows 1:j-1 of current block column
C
            CALL ZTRMM('Left','Upper','No transpose',DIAG,J-1,JB,ONE,A,
     *                 LDA,A(1,J),LDA)
            CALL ZTRSM('Right','Upper','No transpose',DIAG,J-1,JB,-ONE,
     *                 A(J,J),LDA,A(1,J),LDA)
C
C           Compute inverse of current diagonal block
C
            CALL F07TWZ('Upper',DIAG,JB,A(J,J),LDA,INFO)
   40    CONTINUE
      ELSE
C
C        Compute inverse of lower triangular matrix
C
         NN = ((N-1)/NB)*NB + 1
         DO 60 J = NN, 1, -NB
            JB = MIN(NB,N-J+1)
            IF (J+JB.LE.N) THEN
C
C              Compute rows j+jb:n of current block column
C
               CALL ZTRMM('Left','Lower','No transpose',DIAG,N-J-JB+1,
     *                    JB,ONE,A(J+JB,J+JB),LDA,A(J+JB,J),LDA)
               CALL ZTRSM('Right','Lower','No transpose',DIAG,N-J-JB+1,
     *                    JB,-ONE,A(J,J),LDA,A(J+JB,J),LDA)
            END IF
C
C           Compute inverse of current diagonal block
C
            CALL F07TWZ('Lower',DIAG,JB,A(J,J),LDA,INFO)
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of F07TWF (ZTRTRI)
C
      END
