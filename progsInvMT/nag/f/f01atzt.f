      SUBROUTINE F01ATZ(M,A,IA,D,K,L,N,J)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C
C     AUXILIARY ROUTINE CALLED BY F01ATF.
C     INTERCHANGES ELEMENTS 1 TO K OF COLUMNS J AND M,
C     AND ELEMENTS L TO N OF ROWS J AND M.
C     .. Scalar Arguments ..
      INTEGER           IA, J, K, L, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F
      INTEGER           I
C     .. Executable Statements ..
      D(M) = J
      IF (J.EQ.M) GO TO 60
      DO 20 I = 1, K
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   20 CONTINUE
      IF (L.GT.N) GO TO 60
      DO 40 I = L, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
   60 RETURN
      END
