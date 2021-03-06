      SUBROUTINE F04JGT(N,X,SCALE,SUMSQ,TINY,UNDFLW)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (SCLSQS)
C
C     F04JGT RETURNS THE VALUES SCL AND SUM SUCH THAT
C
C     (SCL**2)*SUM = X(1)**2+X(2)**2+...+X(N)**2+(SCALE**2)*SUMSQ .
C
C     SCL IS OVERWRITTEN ON SCALE AND SUM IS OVERWRITTEN ON SUMSQ.
C
C     THE SUPPLIED VALUE OF SUMSQ IS ASSUMED TO BE AT LEAST
C     UNITY IN WHICH CASE SUM WILL SATISFY THE BOUNDS
C
C     1.0 .LE. SUM .LE. SUMSQ+N .
C
C     ONLY ONE PASS THROUGH THE VECTOR X IS MADE.
C
C     TINY MUST BE SUCH THAT
C
C     TINY = SQRT(X02AMF) ,
C
C     WHERE X02AMF IS THE SMALL VALUE RETURNED FROM ROUTINE X02AMF.
C
C     UNDFLW MUST BE THE VALUE RETURNED BY X02DAF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SCALE, SUMSQ, TINY
      INTEGER           N
      LOGICAL           UNDFLW
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSXI, Q
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      IF (UNDFLW) GO TO 60
C
      DO 40 I = 1, N
         IF (X(I).EQ.0.0D0) GO TO 40
C
         ABSXI = ABS(X(I))
         IF (SCALE.GE.ABSXI) GO TO 20
C
         SUMSQ = 1.0D0 + SUMSQ*(SCALE/ABSXI)**2
         SCALE = ABSXI
         GO TO 40
C
   20    SUMSQ = SUMSQ + (ABSXI/SCALE)**2
C
   40 CONTINUE
C
      RETURN
C
   60 DO 100 I = 1, N
         IF (X(I).EQ.0.0D0) GO TO 100
C
         ABSXI = ABS(X(I))
         Q = 0.0D0
         IF (SCALE.LT.ABSXI) GO TO 80
C
         IF (SCALE.GT.TINY) Q = SCALE*TINY
         IF (ABSXI.GE.Q) SUMSQ = SUMSQ + (ABSXI/SCALE)**2
         GO TO 100
C
   80    IF (ABSXI.GT.TINY) Q = ABSXI*TINY
         IF (SCALE.GE.Q) SUMSQ = 1.0D0 + SUMSQ*(SCALE/ABSXI)**2
         SCALE = ABSXI
C
  100 CONTINUE
C
      RETURN
      END
