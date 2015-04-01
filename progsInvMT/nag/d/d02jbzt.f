      SUBROUTINE D02JBZ(X,I,V,A,IA,IA1,RHS,BDYC,BC)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     BC, BDYC
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RHS, X
      INTEGER           I, IA, IA1, V
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,IA1)
C     .. Subroutine Arguments ..
      EXTERNAL          BC, BDYC
C     .. Scalars in Common ..
      DOUBLE PRECISION  X0, X1
      INTEGER           N
C     .. Local Scalars ..
      INTEGER           J
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AD02JB/X0, X1, N
C     .. Executable Statements ..
      CALL BC(I,J,RHS)
      IF (J.EQ.0 .OR. ABS(J).GT.N) RETURN
      X = X1
      IF (J.LT.0) X = X0
      J = ABS(J)
      A(J,1) = 1.D0
      RETURN
      END
