      SUBROUTINE D02HAW(N,X,Y,F,CF,CF1,M,P)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     AUXILIARY ODE ROUTINE
C     CF1, CF
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  F(N), P(M), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          CF, CF1
C     .. Executable Statements ..
      CALL CF(X,Y,F)
      RETURN
      END
