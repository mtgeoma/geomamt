      SUBROUTINE D02JAY(X,I,A,IA,IA1,RHS,COEFF,CF)
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
      INTEGER           M
C     .. Local Scalars ..
      INTEGER           J, M1
C     .. Common blocks ..
      COMMON            /BD02JA/M
C     .. Executable Statements ..
      M1 = M + 1
      DO 20 J = 1, M1
         A(1,J) = CF(J,X)
   20 CONTINUE
      RHS = CF(0,X)
      RETURN
      END
