      DOUBLE PRECISION FUNCTION G01DDZ(C,N,X)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     THIS FUNCTION SUBPROGRAM CALCULATES THE SUM FROM 0 TO N OF
C
C                           I
C                   C(I) * X
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          N
C     .. Array Arguments ..
      DOUBLE PRECISION                 C(0:N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SUM
      INTEGER                          I, J
C     .. Executable Statements ..
      SUM = C(N)
      I = N
      DO 20 J = 1, N
         SUM = SUM*X + C(I-1)
         I = I - 1
   20 CONTINUE
      G01DDZ = SUM
C
      RETURN
      END
