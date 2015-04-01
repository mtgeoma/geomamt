      SUBROUTINE G08AGZ(N,IR,IV,P)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     GIVEN THE RANKS OF THE N DATA OBSERVATIONS IN THE VECTOR, IR,
C     AND THE TEST STATISTIC, IV, G08AGZ FINDS THE LOWER TAIL
C     PROBABILITY, P, USING AN ALGORITHM PRESENTED BY NEUMANN IN
C     THE STATISTICAL SOFTWARE NEWSLETTER, DECEMBER 1988.
C     THE METHOD IS STILL EXACT FOR THE CASE OF TIES.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IV, N
C     .. Array Arguments ..
      INTEGER           IR(N)
C     .. Local Scalars ..
      INTEGER           J, K, LIM, SHIFT, UPPER
C     .. Local Arrays ..
      DOUBLE PRECISION  PROB(3241)
C     .. Executable Statements ..
      DO 20 J = 1, IV + 1
         PROB(J) = 1.0D0
   20 CONTINUE
      UPPER = 0
      DO 60 J = 1, N
         SHIFT = IR(J)
         IF (SHIFT.GE.(IV+1)) THEN
            GO TO 80
         ELSE
            UPPER = UPPER + SHIFT
            LIM = UPPER + 1
            IF (UPPER.GT.IV) LIM = IV + 1
            DO 40 K = LIM, 1, -1
               PROB(K) = 0.5D0*PROB(K)
               IF (SHIFT.LE.K-1) PROB(K) = PROB(K) + 0.5D0*PROB(K-SHIFT)
   40       CONTINUE
         END IF
   60 CONTINUE
      GO TO 100
   80 PROB(IV+1) = PROB(IV+1)/(2.0D0**(N-J+1))
  100 P = PROB(IV+1)
      IF (P.LT.0.0D0) P = 0.0D0
      IF (P.GT.1.0D0) P = 1.0D0
C
      END
