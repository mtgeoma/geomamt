      SUBROUTINE G02AAX(UPLO,N,A,B)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     PLACES THE TRANSPOSE OF A TRIANGULAR MATRIX STORED
C     COLUMN-WISE IN A INTO B
C     UPLO - UPPER ('U') OR LOWER ('L') TRIANGULAR
C     N - SIZE OF MATRIX
C
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER*1       UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N*(N+1)/2), B(N*(N+1)/2)
C     .. Local Scalars ..
      INTEGER           I, IERROR, IJ, J
C     .. Executable Statements ..
      IERROR = 0
      IJ = 0
      IF (UPLO.EQ.'U' .OR. UPLO.EQ.'u') THEN
         DO 40 I = 1, N
            DO 20 J = I, N
               IJ = IJ + 1
               B(IJ) = A((J-1)*J/2+I)
   20       CONTINUE
   40    CONTINUE
      ELSE IF (UPLO.EQ.'L' .OR. UPLO.EQ.'l') THEN
         DO 80 I = 1, N
            DO 60 J = 1, I
               IJ = IJ + 1
               B(IJ) = A((J-1)*N-(J-1)*(J-2)/2+I-J+1)
   60       CONTINUE
   80    CONTINUE
      ELSE
         IERROR = 1
         RETURN
      END IF
      END
