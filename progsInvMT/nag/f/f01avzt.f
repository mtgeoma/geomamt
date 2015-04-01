      SUBROUTINE F01AVZ(M,AR,IAR,AI,IAI,D,K,L,N,J)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C
C     AUXILIARY ROUTINE CALLED BY F01AVF
C     INTERCHANGES ELEMENTS 1 TO K OF COLUMNS J AND M,
C     AND ELEMENTS L TO N OF ROWS J AND M.
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, J, K, L, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), D(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F
      INTEGER           I
C     .. Executable Statements ..
      D(M) = J
      IF (J.EQ.M) GO TO 60
      DO 20 I = 1, K
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   20 CONTINUE
      IF (L.GT.N) GO TO 60
      DO 40 I = L, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
   60 RETURN
      END
