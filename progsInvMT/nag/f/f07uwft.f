      SUBROUTINE F07UWF(UPLO,DIAG,N,AP,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZTPTRI(UPLO,DIAG,N,AP,INFO)
C
C  Purpose
C  =======
C
C  ZTPTRI computes the inverse of a complex upper or lower triangular
C  matrix A stored in packed format.
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
C  AP      (input/output) COMPLEX array, dimension (N*(N+1)/2)
C          On entry, the upper or lower triangular matrix A, stored
C          columnwise in a linear array:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i.
C          If DIAG = 'U', the diagonal elements of A are not referenced
C          and are assumed to be 1.
C          On exit, the (triangular) inverse of the original matrix, in
C          the same packed storage format.
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
      INTEGER           INFO, N
      CHARACTER         DIAG, UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*)
C     .. Local Scalars ..
      COMPLEX*16        AJJ
      INTEGER           J, JC, JCLAST
      LOGICAL           NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          ZSCAL, ZTPMV, F06AAZ
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
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07UWF/ZTPTRI',-INFO)
         RETURN
      END IF
C
C     Check for singularity.
C
      IF (NOUNIT) THEN
         IF (UPPER) THEN
            JC = 1
            DO 20 INFO = 1, N
               IF (AP(JC+INFO-1).EQ.ZERO) RETURN
               JC = JC + INFO
   20       CONTINUE
         ELSE
            JC = 1
            DO 40 INFO = 1, N
               IF (AP(JC).EQ.ZERO) RETURN
               JC = JC + N - INFO + 1
   40       CONTINUE
         END IF
      END IF
      INFO = 0
C
      IF (UPPER) THEN
C
C        Compute inverse of upper triangular matrix.
C
         JC = 1
         DO 60 J = 1, N
            IF (NOUNIT) THEN
               AP(JC+J-1) = ONE/AP(JC+J-1)
               AJJ = -AP(JC+J-1)
            ELSE
               AJJ = -ONE
            END IF
C
C           Compute elements 1:j-1 of j-th column.
C
            CALL ZTPMV('Upper','No transpose',DIAG,J-1,AP,AP(JC),1)
            CALL ZSCAL(J-1,AJJ,AP(JC),1)
            JC = JC + J
   60    CONTINUE
C
      ELSE
C
C        Compute inverse of lower triangular matrix.
C
         JC = N*(N+1)/2
         DO 80 J = N, 1, -1
            IF (NOUNIT) THEN
               AP(JC) = ONE/AP(JC)
               AJJ = -AP(JC)
            ELSE
               AJJ = -ONE
            END IF
            IF (J.LT.N) THEN
C
C              Compute elements j+1:n of j-th column.
C
               CALL ZTPMV('Lower','No transpose',DIAG,N-J,AP(JCLAST),
     *                    AP(JC+1),1)
               CALL ZSCAL(N-J,AJJ,AP(JC+1),1)
            END IF
            JCLAST = JC
            JC = JC - N + J - 2
   80    CONTINUE
      END IF
C
      RETURN
C
C     End of F07UWF (ZTPTRI)
C
      END
