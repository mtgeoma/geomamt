      SUBROUTINE G04EAX(N,IFACT,LEVELS,MEAN,FIRST,X,LDX)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.

C     .. Scalar Arguments ..
      INTEGER           LDX, LEVELS, N
      LOGICAL           FIRST, MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  X(LDX,*)
      INTEGER           IFACT(N)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      IF (MEAN) THEN
         DO 40 J = 1, LEVELS
            DO 20 I = 1, N
               X(I,J) = 0.0D0
   20       CONTINUE
   40    CONTINUE
         DO 60 I = 1, N
            X(I,IFACT(I)) = 1.0D0
   60    CONTINUE
      ELSE IF (FIRST) THEN
         DO 100 J = 1, LEVELS - 1
            DO 80 I = 1, N
               X(I,J) = 0.0D0
   80       CONTINUE
  100    CONTINUE
         DO 120 I = 1, N
            IF (IFACT(I).GT.1) X(I,IFACT(I)-1) = 1.0D0
  120    CONTINUE
      ELSE
         DO 160 J = 1, LEVELS - 1
            DO 140 I = 1, N
               X(I,J) = 0.0D0
  140       CONTINUE
  160    CONTINUE
         DO 180 I = 1, N
            IF (IFACT(I).LT.LEVELS) X(I,IFACT(I)) = 1.0D0
  180    CONTINUE
      END IF
      RETURN
      END
