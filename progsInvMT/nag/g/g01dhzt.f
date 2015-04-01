      SUBROUTINE G01DHZ(SCORES,N,R)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER         SCORES
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIV, THIRD, XN, XN1
      INTEGER           I, IF2, M
C     .. External Functions ..
      DOUBLE PRECISION  G01CEF, G01DAY
      EXTERNAL          G01CEF, G01DAY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      XN = DBLE(N)
      IF (SCORES.EQ.'R') THEN
         DO 20 I = 1, N
            R(I) = DBLE(I)
   20    CONTINUE
      ELSE IF (SCORES.EQ.'B') THEN
         DIV = XN + 0.25D0
         DO 40 I = 1, N
            R(I) = G01CEF((DBLE(I)-0.375D0)/DIV,IF2)
   40    CONTINUE
      ELSE IF (SCORES.EQ.'T') THEN
         THIRD = 1.0D0/3.0D0
         DIV = XN + THIRD
         DO 60 I = 1, N
            R(I) = G01CEF((DBLE(I)-THIRD)/DIV,IF2)
   60    CONTINUE
      ELSE IF (SCORES.EQ.'V') THEN
         DIV = XN + 1.0D0
         DO 80 I = 1, N
            R(I) = G01CEF(DBLE(I)/DIV,IF2)
   80    CONTINUE
      ELSE IF (SCORES.EQ.'N') THEN
         M = N/2
         DO 100 I = 1, M
            R(I) = G01DAY(I,N,IF2)
  100    CONTINUE
         IF (2*M.LT.N) R(M+1) = 0.0D0
         DO 120 I = M + 1, N
            R(I) = -R(N+1-I)
  120    CONTINUE
      ELSE IF (SCORES.EQ.'S') THEN
         R(1) = 1.0D0/XN
         XN1 = XN + 1.0D0
         DO 140 I = 2, N
            R(I) = R(I-1) + 1.0D0/(XN1-DBLE(I))
  140    CONTINUE
      END IF
      RETURN
      END
