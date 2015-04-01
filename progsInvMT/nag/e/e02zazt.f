      INTEGER FUNCTION E02ZAZ(N,X,T)
C     MARK 6 RELEASE. NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     GIVEN REAL T AND A REAL ARRAY X OF DIMENSION
C     AT LEAST N - 4, X(5) TO X(N-4) BEING IN NON-
C     DECREASING ORDER, THIS FUNCTION DETERMINES THE
C     INTERVAL X(I+3) TO X(I+4) WHICH CONTAINS T.
C     SPECIFICALLY, E02ZAZ IS ASSIGNED A VALUE I SUCH
C     THAT X(I+3) .LE. T .LT. X(I+4), UNLESS T .LT. X(5)
C     WHEN E02ZAZ = 1, OR T .GE. X(N-4) WHEN E02ZAZ = N-7.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        T
      INTEGER                 N
C     .. Array Arguments ..
      DOUBLE PRECISION        X(N)
C     .. Local Scalars ..
      INTEGER                 I, L, U
C     .. Executable Statements ..
      L = 4
      U = N - 3
   20 I = (L+U)/2
      IF (U-L.LE.1) GO TO 60
      IF (T.GE.X(I)) GO TO 40
      U = I
      GO TO 20
   40 L = I
      GO TO 20
   60 E02ZAZ = U - 4
      RETURN
      END
