      SUBROUTINE G02AAY(UPLO,DIAG,N,A)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     PRE-MULTIPLIES A TRIANGULAR MATRIX (A) STORED COLUMN-WISE
C     BY ITS TRANSPOSE
C     UPLO - UPPER ('U') OR LOWER ('L') TRIANGULAR MATRIX
C     DIAG - NON-UNIT ('N') OR UNIT ('U') TRIANGULAR
C     N - SIZE OF MATRIX
C     WK - WORKSPACE
C
C
C     CHECK FOR ERRORS IN INPUT
C
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER*1       DIAG, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N*(N+1)/2)
C     .. Local Scalars ..
      INTEGER           I, IERROR, J, K
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DTPMV
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
C
C     CALCULATE A'A
C
         IF (UPLO.EQ.'U' .OR. UPLO.EQ.'u') THEN
            K = (N-1)*N/2 + 1
            IF (DIAG.EQ.'N' .OR. DIAG.EQ.'n') THEN
               DO 20 I = N - 1, 1, -1
                  J = K + I
                  A(J) = A(J)*A(J) + DDOT(I,A(K),1,A(K),1)
                  CALL DTPMV('U','T','N',I,A,A(K),1)
                  K = K - I
   20          CONTINUE
               A(1) = A(1)*A(1)
            ELSE
               DO 40 I = N - 1, 1, -1
                  J = K + I
                  A(J) = DDOT(I,A(K),1,A(K),1) + 1.0D0
                  CALL DTPMV('U','T','U',I,A,A(K),1)
                  K = K - I
   40          CONTINUE
               A(1) = 1.0D0
            END IF
         ELSE
            J = 1
            IF (DIAG.EQ.'N' .OR. DIAG.EQ.'n') THEN
               DO 60 I = N - 1, 1, -1
                  K = J + 1
                  A(J) = A(J)*A(J) + DDOT(I,A(K),1,A(K),1)
                  J = K + I
                  CALL DTPMV('L','T','N',I,A(J),A(K),1)
   60          CONTINUE
               A(J) = A(J)*A(J)
            ELSE
               DO 80 I = N - 1, 1, -1
                  K = J + 1
                  A(J) = DDOT(I,A(K),1,A(K),1) + 1.0D0
                  J = K + I
                  CALL DTPMV('L','T','U',I,A(J),A(K),1)
   80          CONTINUE
               A(J) = 1.0D0
            END IF
         END IF
      END IF
      END
