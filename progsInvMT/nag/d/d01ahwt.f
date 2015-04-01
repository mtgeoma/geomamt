      SUBROUTINE D01AHW(F,V,A,B,R1,VAL,IT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     CALCULATES THE VALUE OF THE TRANSFORMED INTEGRAND.
C     A RANDOM TRANSFORMATION (IR NON-ZERO, AND CONTROLLED BY ALP)
C     AND POSSIBLY A SINGULARITY WEAKENING TRANSFORMATION
C     (IT NON-ZERO, AND CONTROLLED BY NT) WILL BE APPLIED.
C     WITH A POSSIBLE SINGULARITY AT ENDPOINT  R1  THE
C     TRANSFORMATION OF VARIABLE IS
C
C                X = (T-A)**(NT+1)/(B-A)**NT+A,  R1=A
C          OR    X = (T-B)**(NT+1)/(A-B)**NT+B,  R1=B
C
C     THE RANDOM TRANSFORMATION IS -
C
C                X = A*B*ALP+(1-ALP*(A+B))*T+ALP*T**2
C
C     WHERE  ALP*(B-A)  IS    RANDOM IN  (-.01,.01).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, R1, V, VAL
      INTEGER           IT
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  AFLOW, ALP, AV, EPMACH, UFLOW
      INTEGER           IR, NT
C     .. Local Scalars ..
      DOUBLE PRECISION  ELIM, GM, PART, Q, RK, S, SLOPE, SLOPER, T, W, X
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, DBLE
C     .. Common blocks ..
      COMMON            /CD01AH/ALP, AV, NT, IR
      COMMON            /DD01AH/EPMACH, UFLOW, AFLOW
C     .. Executable Statements ..
      SLOPER = 1.0D0
      T = V
      IF (IR.EQ.0) GO TO 20
C
C     RANDOM TRANSFORMATION
      PART = 1.0D0 - ALP*(A+B)
      T = ALP*A*B + (PART+ALP*V)*V
      SLOPER = PART + 2.0D0*ALP*V
   20 IF (IT.EQ.0) GO TO 100
      IF (B.EQ.R1) GO TO 80
C
C     LEFT ENDPOINT PEAK
      RK = B - A
      W = A
   40 S = T - W
      GM = 0.0D0
      ELIM = AFLOW/DBLE(NT)
      Q = S/RK
      IF (ABS(Q).GT.EXP(ELIM)) GM = Q**NT
      X = 0.0D0
      IF (S.EQ.0.0D0) GO TO 60
      IF (ABS(GM).GT.UFLOW/ABS(S)) X = S*GM
   60 SLOPE = DBLE(NT+1)*GM*SLOPER
      GO TO 120
C
C     RIGHT ENDPOINT PEAK
   80 RK = A - B
      W = B
      GO TO 40
  100 X = T
      W = 0.0D0
      SLOPE = SLOPER
  120 VAL = SLOPE*F(X+W)
      RETURN
      END
