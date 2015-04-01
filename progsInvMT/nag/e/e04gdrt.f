      SUBROUTINE E04GDR(M,N,A,LA,X,Y)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     THIS ROUTINE CALCULATES A*X + Y AND OVERWRITES THE RESULT
C     ON Y.
C
C     SVEN HAMMARLING
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           LA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA,N), X(N), Y(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  XJ
      INTEGER           I, J
C     .. Executable Statements ..
      DO 40 J = 1, N
         XJ = X(J)
         DO 20 I = 1, M
            Y(I) = Y(I) + XJ*A(I,J)
   20    CONTINUE
   40 CONTINUE
      RETURN
C
C     END OF E04GDR   (AXPLSY)
C
      END
