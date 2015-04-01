      SUBROUTINE D02SAY(N,X,Y,F,CF,CF1,M,P)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11 REVISED. IER-419 (FEB 1984).
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
C     .. Scalars in Common ..
      INTEGER           I, ICASE, IW2
C     .. Arrays in Common ..
      DOUBLE PRECISION  COUT(2), W(7)
      INTEGER           IW1(4)
C     .. Common blocks ..
      COMMON            /AD02SA/W, IW1, IW2, ICASE
      COMMON            /BD02SA/COUT, I
C     .. Executable Statements ..
      IF (ICASE.GT.1) GO TO 20
      CALL CF1(X,Y,F,N,P,M,I)
      RETURN
   20 CALL CF(X,Y,F,P)
      RETURN
      END
