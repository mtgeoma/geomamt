      SUBROUTINE D02TGW(X,I,V,A,IA,IA1,RHS,BDYC,BC)
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
C     .. Executable Statements ..
      CALL BDYC(X,I,V,A,IA,IA1,RHS)
      RETURN
      END
