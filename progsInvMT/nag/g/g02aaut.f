      SUBROUTINE G02AAU(I,J,N,A)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     SWOPS ROWS AND COLUMNS OF A REAL SYMMETRIC MATRIX
C     STORED IN PACKED UPPER TRIANGULAR FORM BY COLUMNS.
C
C     .. Scalar Arguments ..
      INTEGER           I, J, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N*(N+1)/2)
C     .. Local Scalars ..
      DOUBLE PRECISION  SWAP
      INTEGER           K, L, LK, U, UK
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
      IF (N.GT.1) THEN
         IF (I.NE.J) THEN
            U = MAX(I,J)
            L = MIN(I,J)
            LK = (L-1)*L/2
            UK = (U-1)*U/2
            DO 20 K = 1, L - 1
               LK = LK + 1
               UK = UK + 1
               SWAP = A(LK)
               A(LK) = A(UK)
               A(UK) = SWAP
   20       CONTINUE
            LK = LK + 1
            UK = UK + 1
            SWAP = A(LK)
            A(LK) = A(UK+U-L)
            A(UK+U-L) = SWAP
            DO 40 K = L, U - 2
               LK = LK + K
               UK = UK + 1
               SWAP = A(LK)
               A(LK) = A(UK)
               A(UK) = SWAP
   40       CONTINUE
            LK = LK + U - 1
            UK = UK + 1
            DO 60 K = U, N - 1
               LK = LK + K
               UK = UK + K
               SWAP = A(LK)
               A(LK) = A(UK)
               A(UK) = SWAP
   60       CONTINUE
         END IF
      END IF
      END
