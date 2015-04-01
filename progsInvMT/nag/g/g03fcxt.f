      SUBROUTINE G03FCX(N,M,X,D,NN,SRD)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Calculate Euclidean distance matrix for X
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NN
      LOGICAL           SRD
C     .. Array Arguments ..
      DOUBLE PRECISION  D(NN), X(N,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  S
      INTEGER           I, II, J, K
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      II = 0
      DO 60 I = 2, N
         DO 40 J = 1, I - 1
            II = II + 1
            S = 0.0D0
            DO 20 K = 1, M
               S = S + (X(I,K)-X(J,K))**2
   20       CONTINUE
            IF (SRD) S = SQRT(S)
            D(II) = S
   40    CONTINUE
   60 CONTINUE
      RETURN
      END
