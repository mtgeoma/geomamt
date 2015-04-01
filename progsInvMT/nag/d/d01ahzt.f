      DOUBLE PRECISION FUNCTION D01AHZ(S1,S3,AL,AR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBDIVIDE IN RATIO 1/2 IF INTEGRAND IS STEEPER ON LEFT OR
C                        2/1 IF STEEPER ON RIGHT
C                        1/1 OTHERWISE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 AL, AR, S1, S3
C     .. Executable Statements ..
      IF (AL-AR) 20, 40, 60
   20 D01AHZ = (S1+2.0D0*S3)/3.0D0
      RETURN
   40 D01AHZ = (S1+S3)/2.0D0
      RETURN
   60 D01AHZ = (2.0D0*S1+S3)/3.0D0
      RETURN
      END
