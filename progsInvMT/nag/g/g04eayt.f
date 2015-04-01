      SUBROUTINE G04EAY(N,IFACT,LEVELS,REP,X,LDX)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.

C     .. Scalar Arguments ..
      INTEGER           LDX, LEVELS, N
C     .. Array Arguments ..
      DOUBLE PRECISION  REP(LEVELS), X(LDX,*)
      INTEGER           IFACT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, TEMP
      INTEGER           I, J, K
C     .. Executable Statements ..
      SUM = REP(1)
      DO 20 I = 2, LEVELS
         TEMP = REP(I)
         REP(I) = SUM/TEMP
         SUM = SUM + TEMP
   20 CONTINUE
      DO 80 I = 1, N
         K = IFACT(I)
         DO 40 J = 1, K - 1
            X(I,J) = 0.0D0
   40    CONTINUE
         IF (K.GT.1) X(I,K-1) = REP(K)
         DO 60 J = K, LEVELS - 1
            X(I,J) = -1.0D0
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
