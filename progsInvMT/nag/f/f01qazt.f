      SUBROUTINE F01QAZ(N,Z,Z1,X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (HOUSH1)
C
C     F01QAZ RETURNS THE N ELEMENT VECTOR
C
C     Y = (I-(1/Z(1))*Z*(Z**T))*X ,
C
C     WHERE X AND Z ARE N ELEMENT VECTORS.
C
C     Y IS OVERWRITTEN ON X.
C
C     THE VALUE OF Z(1) MUST ACTUALLY BE SUPPLIED IN Z1.
C     THE ELEMENT Z(1) IS NOT REFERENCED.
C
C
C     N MUST BE AT LEAST 1. IF N=1 THEN AN IMMEDIATE RETURN TO
C     THE CALLING PROGRAM IS MADE.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  Z1
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Z(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D
      LOGICAL           UNDFLW
C     .. External Functions ..
      DOUBLE PRECISION  F01QAX
      LOGICAL           X02DAF
      EXTERNAL          F01QAX, X02DAF
C     .. External Subroutines ..
      EXTERNAL          F01QAW
C     .. Executable Statements ..
      IF (N.EQ.1) RETURN
C
      UNDFLW = X02DAF(0.0D0)
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     THE CALL TO F01QAX CAN BE REPLACED BY THE FOLLOWING IN-LINE
C     CODE, PROVIDED THAT NO PRECAUTIONS AGAINST UNDERFLOW
C     ARE REQUIRED
C
C     D = X(1)*Z1
C     DO 20 I=2,N
C        D = D + Z(I)*X(I)
C     20 CONTINUE
C
C     IN THIS CASE THE DECLARATION
C
C     REAL F01QAX
C
C     MUST ALSO BE REMOVED.
C
      D = F01QAX(N-1,N-1,X(1)*Z1,.TRUE.,Z(2),X(2),UNDFLW)
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      X(1) = X(1) - D
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     THE CALL TO F01QAW CAN BE REPLACED BY THE FOLLOWING IN-LINE
C     CODE, PROVIDED THAT NO PRECAUTIONS AGAINST UNDERFLOW
C     ARE REQUIRED
C
C     D = D/Z1
C     DO 40 I=2,N
C        X(I) = X(I) - D*Z(I)
C     40 CONTINUE
C
      CALL F01QAW(N-1,X(2),D/Z1,Z(2),UNDFLW)
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      RETURN
      END
