      DOUBLE PRECISION FUNCTION G13BCY(A,N)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUMS THE ELEMENTS OF ARRAY A
C
C     .. Scalar Arguments ..
      INTEGER                          N
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(N)
C     .. Local Scalars ..
      INTEGER                          I
C     .. Executable Statements ..
      G13BCY = 0.0D0
      DO 20 I = 1, N
         G13BCY = G13BCY + A(I)
   20 CONTINUE
      RETURN
      END
