      SUBROUTINE F01QAW(N,X,XMUL,Y,UNDFLW)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (ROWOP2)
C
C     F01QAW RETURNS THE N ELEMENT VECTOR Z GIVEN BY
C
C     Z = X - XMUL*Y ,
C
C     WHERE XMUL IS A REAL VALUE AND X AND Y ARE N ELEMENT
C     VECTORS.
C
C     Z IS OVERWRITTEN ON X.
C
C     N MUST BE AT LEAST 1.
C
C     UNDFLW MUST BE THE VALUE RETURNED BY X02DAF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMUL
      INTEGER           N
      LOGICAL           UNDFLW
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMUL, W
      INTEGER           I
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      IF (XMUL.EQ.0.0D0) RETURN
C
      IF (UNDFLW) GO TO 60
C
   20 DO 40 I = 1, N
         X(I) = X(I) - XMUL*Y(I)
   40 CONTINUE
C
      RETURN
C
   60 AMUL = ABS(XMUL)
      IF (AMUL.GE.1.0D0) GO TO 20
      W = X02AMF()/AMUL
C
      DO 80 I = 1, N
         IF (ABS(Y(I)).LT.W) GO TO 80
         X(I) = X(I) - XMUL*Y(I)
   80 CONTINUE
C
      RETURN
      END
