      DOUBLE PRECISION FUNCTION F01QAX(NR,N,V,PLUS,X,Y,UNDFLW)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (DOTPRD)
C
C     F01QAX RETURNS THE VALUE
C
C     F01QAX = ( V + (X**T)*Y , WHEN PLUS = .TRUE.
C              (
C              ( V - (X**T)*Y , WHEN PLUS = .FALSE. ,
C
C     WHERE V IS A REAL VALUE AND X AND Y ARE N ELEMENT VECTORS.
C
C     IF N IS LESS THAN UNITY THEN F01QAX IS RETURNED AS V.
C
C     NR MUST BE AT LEAST MAX(1,N).
C
C     UNDFLW MUST BE THE VALUE RETURNED BY X02DAF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 V
      INTEGER                          N, NR
      LOGICAL                          PLUS, UNDFLW
C     .. Array Arguments ..
      DOUBLE PRECISION                 X(NR), Y(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSXI, SMALL, SUM
      INTEGER                          I
C     .. External Functions ..
      DOUBLE PRECISION                 X02AMF
      EXTERNAL                         X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Executable Statements ..
      SUM = V
C
      IF (N.LT.1) GO TO 80
C
      IF (UNDFLW) GO TO 100
C
      IF (PLUS) GO TO 40
C
      DO 20 I = 1, N
         SUM = SUM - X(I)*Y(I)
   20 CONTINUE
C
      F01QAX = SUM
C
      RETURN
C
   40 DO 60 I = 1, N
         SUM = SUM + X(I)*Y(I)
   60 CONTINUE
C
   80 F01QAX = SUM
C
      RETURN
C
  100 SMALL = X02AMF()
C
      IF (PLUS) GO TO 160
C
      DO 140 I = 1, N
         ABSXI = ABS(X(I))
         IF (ABSXI.GE.1.0D0 .OR. ABSXI.EQ.0.0D0) GO TO 120
         IF (ABS(Y(I)).LT.SMALL/ABSXI) GO TO 140
  120    SUM = SUM - X(I)*Y(I)
  140 CONTINUE
C
      F01QAX = SUM
C
      RETURN
C
  160 DO 200 I = 1, N
         ABSXI = ABS(X(I))
         IF (ABSXI.GE.1.0D0 .OR. ABSXI.EQ.0.0D0) GO TO 180
         IF (ABS(Y(I)).LT.SMALL/ABSXI) GO TO 200
  180    SUM = SUM + X(I)*Y(I)
  200 CONTINUE
C
      F01QAX = SUM
C
      RETURN
      END
