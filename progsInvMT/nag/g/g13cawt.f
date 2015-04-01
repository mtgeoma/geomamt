      DOUBLE PRECISION FUNCTION G13CAW(IW,S3,PI)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13CAW CALCULATES THE VALUE OF THE LAG WINDOW
C
C     IW    - SHAPE OF LAG WINDOW
C     1 RECTANGULAR   2 BARTLETT
C     3 TUKEY         4 PARZEN
C     S3    - ARGUMENT OF WINDOW FUNCTION
C     PI    - CONSTANT PI
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 PI, S3
      INTEGER                          IW
C     .. Local Scalars ..
      DOUBLE PRECISION                 S2
C     .. Intrinsic Functions ..
      INTRINSIC                        COS
C     .. Executable Statements ..
      GO TO (20,40,60,80) IW
   20 S2 = 1.0D0
      GO TO 120
   40 S2 = 1.0D0 - S3
      GO TO 120
   60 S2 = (1.0D0+COS(PI*S3))/2.0D0
      GO TO 120
   80 IF (S3.GT.0.5D0) GO TO 100
      S2 = 1.0D0 - 6.0D0*S3*S3*(1.0D0-S3)
      GO TO 120
  100 S2 = 2.0D0*(1.0D0-S3)**3
  120 G13CAW = S2
      RETURN
      END
