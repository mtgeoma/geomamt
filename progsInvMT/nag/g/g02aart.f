      DOUBLE PRECISION FUNCTION G02AAR(UPLO,N,A)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes the trace of a packed (symmetric) matrix
C
C     .. Scalar Arguments ..
      INTEGER                          N
      CHARACTER                        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SUM
      INTEGER                          I, II
C     .. Executable Statements ..
      IF (N.GT.0) THEN
         IF (N.EQ.1) THEN
            G02AAR = A(1)
         ELSE IF (UPLO.EQ.'U' .OR. UPLO.EQ.'u') THEN
            SUM = A(1)
            II = 1
            DO 20 I = 2, N
               II = II + I
               SUM = SUM + A(II)
   20       CONTINUE
            G02AAR = SUM
         ELSE IF (UPLO.EQ.'L' .OR. UPLO.EQ.'l') THEN
            SUM = A(1)
            II = 1
            DO 40 I = 1, N - 1
               II = II + N - I + 1
               SUM = SUM + A(II)
   40       CONTINUE
            G02AAR = SUM
         END IF
      END IF
      RETURN
      END
