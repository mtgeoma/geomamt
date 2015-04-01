      SUBROUTINE G13BCZ(X,Y,NXY,NL,S,R0,R,STAT,IERROR)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     G13BCZ CARRIES OUT THE CALCULATIONS FOR G13BCF
C
C     X     - X SERIES ARRAY
C     Y     - Y SERIES ARRAY
C     NXY   - LENGTH OF X AND Y SERIES
C     NL    - NO OF LAGS
C     S     - STANDARD DEVIATION RATIO
C     R0    - LAG ZERO CORRELATION BETWEEN X AND Y SERIES
C     R     - CROSS CORRELATION ARRAY
C     STAT  - TEST STATISTIC FOR LACK OF CROSS CORRELATION
C     IERROR- ERROR INDICATOR
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  R0, S, STAT
      INTEGER           IERROR, NL, NXY
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NL), X(NXY), Y(NXY)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, EPS, SX, SY, XM, YM, Z
      INTEGER           I, J, M
C     .. External Functions ..
      DOUBLE PRECISION  G13BCX, G13BCY, X02AJF
      EXTERNAL          G13BCX, G13BCY, X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      EPS = X02AJF()
      Z = DBLE(NXY)
      XM = G13BCY(X,NXY)/Z
      YM = G13BCY(Y,NXY)/Z
      SX = SQRT(G13BCX(X,X,NXY,XM,XM)/Z)
      IF (SX.LT.EPS) GO TO 60
      SY = SQRT(G13BCX(Y,Y,NXY,YM,YM)/Z)
      IF (SY.LT.EPS) GO TO 60
      C = Z*SX*SY
      IF (C.LT.EPS) GO TO 60
      R0 = G13BCX(X,Y,NXY,XM,YM)/C
      M = NXY
      J = 1
      DO 20 I = 1, NL
         M = M - 1
         J = J + 1
         R(I) = G13BCX(X(1),Y(J),M,XM,YM)/C
   20 CONTINUE
      STAT = 0
      DO 40 I = 1, NL
         STAT = STAT + R(I)*R(I)
   40 CONTINUE
      STAT = STAT*Z
      S = SY/SX
C     SUCCESSFUL EXIT
      RETURN
C     UNSUCCESSFUL EXIT
C     ZERO VARIANCE
   60 IERROR = 2
      RETURN
      END
