      SUBROUTINE G13DLZ(I,K,N,ZOLD,IK,TRAN,Z,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1690 (JUN 1995).
C     .. Scalar Arguments ..
      INTEGER           I, IFAULT, IK, K, N
      CHARACTER*1       TRAN
C     .. Array Arguments ..
      DOUBLE PRECISION  Z(K,N), ZOLD(IK,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  TS
      INTEGER           L
C     .. External Subroutines ..
      EXTERNAL          DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, SQRT
C     .. Executable Statements ..
C
C     Transform each series (if necessary)
C
      IFAULT = 0
      IF (TRAN.EQ.'L' .OR. TRAN.EQ.'l') THEN
C
C        Log transformation
C
         DO 20 L = 1, N
            TS = ZOLD(I,L)
            IF (TS.GT.0.0D0) THEN
               Z(I,L) = LOG(TS)
            ELSE
               IFAULT = 1
               RETURN
            END IF
   20    CONTINUE
      ELSE IF (TRAN.EQ.'S' .OR. TRAN.EQ.'s') THEN
C
C        Square root transformation
C
         DO 40 L = 1, N
            TS = ZOLD(I,L)
            IF (TS.GE.0.0D0) THEN
               Z(I,L) = SQRT(TS)
            ELSE
               IFAULT = 1
               RETURN
            END IF
   40    CONTINUE
      ELSE
C
C        No transformation
C
         CALL DCOPY(N,ZOLD(I,1),IK,Z(I,1),K)
      END IF
C
      RETURN
      END
