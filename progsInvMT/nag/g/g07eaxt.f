      DOUBLE PRECISION FUNCTION G07EAX(LOWER,X,N,XME,WRK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes the Wilcoxon signed rank statistic.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 XME
      INTEGER                          N
      LOGICAL                          LOWER
C     .. Array Arguments ..
      DOUBLE PRECISION                 WRK(2*N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 RS, XF
      INTEGER                          I
C     .. External Subroutines ..
      EXTERNAL                         G08AEZ
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Executable Statements ..
      DO 20 I = 1, N
         WRK(I) = ABS(X(I)-XME)
   20 CONTINUE
      CALL G08AEZ(WRK(1),WRK(N+1),N,1,XF)
      RS = 0.0D0
      IF (LOWER) THEN
         DO 40 I = 1, N
            IF ((X(I)-XME).GE.0.0D0) RS = RS + WRK(N+I)
   40    CONTINUE
      ELSE
         DO 60 I = 1, N
            IF ((X(I)-XME).GT.0.0D0) RS = RS + WRK(N+I)
   60    CONTINUE
      END IF
      G07EAX = RS
C
      RETURN
      END
