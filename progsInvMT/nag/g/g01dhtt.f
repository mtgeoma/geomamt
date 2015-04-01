      DOUBLE PRECISION FUNCTION G01DHT(SCORES,I,N)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER                          I, N
      CHARACTER                        SCORES
C     .. Local Scalars ..
      DOUBLE PRECISION                 THIRD, W, XN
      INTEGER                          IF2
C     .. External Functions ..
      DOUBLE PRECISION                 G01CEF
      EXTERNAL                         G01CEF
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE
C     .. Executable Statements ..
C
      XN = DBLE(N)
      IF (SCORES.EQ.'R') THEN
         G01DHT = DBLE(I)
      ELSE
         IF (SCORES.EQ.'B') THEN
            W = (DBLE(I)-0.375D0)/(XN+0.25D0)
         ELSE IF (SCORES.EQ.'T') THEN
            THIRD = 1.0D0/3.0D0
            W = (DBLE(I)-THIRD)/(XN+THIRD)
         ELSE IF (SCORES.EQ.'V') THEN
            W = DBLE(I)/(XN+1.0D0)
         END IF
         G01DHT = G01CEF(W,IF2)
      END IF
      RETURN
      END
