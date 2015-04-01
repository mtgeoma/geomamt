      SUBROUTINE D02TGV(X,I,A,IA,IA1,RHS,COEFF,CF)
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
C     .. Executable Statements ..
      CALL COEFF(X,I,A,IA,IA1,RHS)
      RETURN
      END
