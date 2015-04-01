      SUBROUTINE C06FQQ(A,M,N)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:M-1,0:N-1)
C     .. Local Scalars ..
      INTEGER           L
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 L = 0, M - 1
         A(L,0) = 0.5D0*A(L,0)
   20 CONTINUE
      IF (MOD(N,2).EQ.0) THEN
         DO 40 L = 0, M - 1
            A(L,N/2) = 0.5D0*A(L,N/2)
   40    CONTINUE
      END IF
C
      RETURN
      END
