      SUBROUTINE G13BET(A,IDA,N,D)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BET TAKES ARRAY OUTPUT BY F04ASF AND
C     CALCULATES ITS DETERMINANT. IT MAKES USE OF THE FACT
C     THAT THE ORIGINAL ARRAY WAS DECOMPOSED INTO A FORM
C     CONTAINING A LOWER TRIANGULAR MATRIX L.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D
      INTEGER           IDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IDA,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  V
      INTEGER           I, IM, J
C     .. Executable Statements ..
      D = A(1,1)
      IF (N.LE.1) GO TO 60
      DO 40 I = 2, N
         V = A(I,I)
         IM = I - 1
         DO 20 J = 1, IM
            V = V - A(I,J)**2
   20    CONTINUE
         D = D*V
   40 CONTINUE
   60 RETURN
      END
