      SUBROUTINE F04AGZ(A,IA,N,P,B)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMPUTES  B = L**(-1) * B  WHERE
C     A HOLDS THE SUBDIAGONAL ELEMENTS OF L AND
C     P HOLDS THE RECIPROCALS OF THE DIAGONAL ELEMENTS OF L.
C
C
C     .. Scalar Arguments ..
      INTEGER           IA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(N), P(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X
      INTEGER           I, J, JP1, NM1
C     .. Executable Statements ..
      NM1 = N - 1
      IF (NM1.EQ.0) GO TO 60
      DO 40 J = 1, NM1
         B(J) = B(J)*P(J)
         JP1 = J + 1
         X = B(J)
         DO 20 I = JP1, N
            B(I) = B(I) - A(I,J)*X
   20    CONTINUE
   40 CONTINUE
   60 B(N) = B(N)*P(N)
      RETURN
      END
