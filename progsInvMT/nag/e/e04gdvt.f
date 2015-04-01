      SUBROUTINE E04GDV(N,A,LA)
C     MARK 11 RE-ISSUE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     RETURNS THE TRANSPOSE OF A
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           LA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  Z
      INTEGER           I, J, JM1
C     .. Executable Statements ..
      IF (N.EQ.1) RETURN
      DO 40 J = 2, N
         JM1 = J - 1
         DO 20 I = 1, JM1
            Z = A(J,I)
            A(J,I) = A(I,J)
            A(I,J) = Z
   20    CONTINUE
   40 CONTINUE
      RETURN
C
C     END OF E04GDV (TRNSPS)
C
      END
