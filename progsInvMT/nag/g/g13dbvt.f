      SUBROUTINE G13DBV(A,IA,B,IB,N)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        COPIES UPPER TRIANGLE OF TOP LEFT (N,N) SQUARE
C        OF MATRIX A INTO FULL SYMMETRIC FORM IN TOP LEFT
C        (N,N) SQUARE IN MATRIX B. A AND B MAY BE
C        IDENTICAL
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      DO 40 I = 1, N
         DO 20 J = I, N
            B(J,I) = A(I,J)
            B(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
      CONTINUE
      RETURN
      END
