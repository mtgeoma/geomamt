      SUBROUTINE G02AAZ(UPLO,DIAG,N,A)
C
C     FINDS THE INVERSE OF A TRIANGULAR MATRIX STORED
C     PACKED COLUMN-WISE IN ARRAY A
C     N - SIZE OF MATRIX
C     UPLO - UPPER ('U') OR LOWER ('L') TRIANGULAR
C     DIAG - NON-UNIT TRIANGULAR ('N') OR UNIT TRIANGULAR ('U')
C     WK - WORKSPACE
C
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER*1       DIAG, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N*(N+1)/2)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, IERROR, J, K
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DTPSV
C     .. Executable Statements ..
      IF (UPLO.NE.'L' .AND. UPLO.NE.'l' .AND. UPLO.NE.'U' .AND. UPLO.NE.
     *    'u') THEN
         IERROR = 1
      ELSE IF (DIAG.NE.'U' .AND. DIAG.NE.'u' .AND. DIAG.NE.'N' .AND.
     *         DIAG.NE.'n') THEN
         IERROR = 2
      ELSE IF (N.LT.0) THEN
         IERROR = 3
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.0) THEN
         RETURN
      ELSE IF (N.GT.0) THEN
         IF (DIAG.EQ.'N' .OR. DIAG.EQ.'n') THEN
            IF (UPLO.EQ.'L' .OR. UPLO.EQ.'l') THEN
               J = 1
               DO 20 I = N - 1, 1, -1
                  TEMP = 1.0D0/A(J)
                  A(J) = TEMP
                  K = J + 1
                  CALL DSCAL(I,-TEMP,A(K),1)
                  J = K + I
                  CALL DTPSV('L','N','N',I,A(J),A(K),1)
   20          CONTINUE
               A(J) = 1.0D0/A(J)
            ELSE
               K = (N-1)*N/2 + 1
               DO 40 I = N - 1, 1, -1
                  J = K + I
                  TEMP = 1.0D0/A(J)
                  A(J) = TEMP
                  CALL DSCAL(I,-TEMP,A(K),1)
                  CALL DTPSV('U','N','N',I,A(1),A(K),1)
                  K = K - I
   40          CONTINUE
               A(1) = 1.0D0/A(1)
            END IF
         ELSE IF (UPLO.EQ.'L' .OR. UPLO.EQ.'l') THEN
            J = 1
            DO 60 I = N - 1, 1, -1
               A(J) = 1.0D0
               K = J + 1
               CALL DSCAL(I,-1.0D0,A(K),1)
               J = K + I
               CALL DTPSV('L','N','U',I,A(J),A(K),1)
   60       CONTINUE
            A(J) = 1.0D0
         ELSE
            K = (N-1)*N/2 + 1
            DO 80 I = N - 1, 1, -1
               J = K + I
               A(J) = 1.0D0
               CALL DSCAL(I,-1.0D0,A(K),1)
               CALL DTPSV('U','N','U',I,A(1),A(K),1)
               K = K - I
   80       CONTINUE
            A(J) = 1.0D0
         END IF
      END IF
      END
