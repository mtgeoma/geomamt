      SUBROUTINE G13BBY(A,NA,B,NB,K,IA,IB)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        ROUTINE TO TRANSFER DATA SLICE FROM B TO A
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, K, NA, NB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NA), B(NB)
C     .. Local Scalars ..
      INTEGER           I
C     .. Executable Statements ..
      IF (K.EQ.0) GO TO 40
      DO 20 I = 1, K
         IA = IA + 1
         IB = IB + 1
         A(IA) = B(IB)
   20 CONTINUE
   40 CONTINUE
      RETURN
      END
