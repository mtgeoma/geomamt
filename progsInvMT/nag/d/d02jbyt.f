      SUBROUTINE D02JBY(X,I,A,IA,IA1,RHS,COEFF,CF)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COEFF
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RHS, X
      INTEGER           I, IA, IA1
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,IA1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF
      EXTERNAL          CF
C     .. Subroutine Arguments ..
      EXTERNAL          COEFF
C     .. Scalars in Common ..
      DOUBLE PRECISION  X0, X1
      INTEGER           N
C     .. Local Scalars ..
      INTEGER           J
C     .. Common blocks ..
      COMMON            /AD02JB/X0, X1, N
C     .. Executable Statements ..
      DO 20 J = 1, N
         A(J,1) = -CF(I,J,X)
   20 CONTINUE
      A(I,2) = 1.D0
      RHS = CF(I,0,X)
      RETURN
      END
