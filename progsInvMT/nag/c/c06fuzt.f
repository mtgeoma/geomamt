      SUBROUTINE C06FUZ(X,Y,XX,YY,M,N)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     C06FUZ transposes a two-dimensional complex array.
C
C     .. Scalar Arguments ..
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(M,N), XX(N,M), Y(M,N), YY(N,M)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      IF (M.GE.N) THEN
         DO 40 J = 1, N
            DO 20 I = 1, M
               XX(J,I) = X(I,J)
               YY(J,I) = Y(I,J)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 I = 1, M
            DO 60 J = 1, N
               XX(J,I) = X(I,J)
               YY(J,I) = Y(I,J)
   60       CONTINUE
   80    CONTINUE
      END IF
      RETURN
      END
