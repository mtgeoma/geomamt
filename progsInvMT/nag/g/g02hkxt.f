      SUBROUTINE G02HKX(EPS,C)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES HUBER'S CONSTANT C
C
C     .. Parameters ..
      INTEGER           MAX
      PARAMETER         (MAX=25)
      DOUBLE PRECISION  EMIN, CMAX, TOL
      PARAMETER         (EMIN=3.5D-6,CMAX=4.0D0,TOL=0.0005D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, EPS
C     .. Local Scalars ..
      DOUBLE PRECISION  ANC, CC, CONST, D, EX, FP, PH, PI, XMAX
      INTEGER           I, IFAULT
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X01AAF, X02AMF
      EXTERNAL          S15ABF, X01AAF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, SQRT
C     .. Executable Statements ..
C
C    CHECK FOR MINIMUM EPS
C
      IF (EPS.LT.EMIN) THEN
         C = CMAX
      ELSE
C
C     SET CONSTANTS
C
         PI = X01AAF(C)
         ANC = 1.0D0/SQRT(2.0D0*PI)
         XMAX = -LOG(X02AMF())
         CONST = (EPS-2.D0)/(1.D0-EPS)/2.D0
C
C     USE NEWTON'S ALGORITHM
C
         C = 1.0D0
         DO 20 I = 1, MAX
            CC = 0.5D0*C*C
            IF (CC.GE.XMAX) THEN
               EX = 0.0D0
            ELSE
               EX = ANC*EXP(-CC)
            END IF
            IFAULT = 1
            PH = S15ABF(C,IFAULT)
            FP = PH + CONST
            D = -EX/FP
            IF (ABS((D-C)/(D+1.0D0)).LT.TOL) GO TO 40
            C = D
   20    CONTINUE
   40    CONTINUE
      END IF
      RETURN
      END
