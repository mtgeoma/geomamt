      SUBROUTINE E02BEV(T,N,K,X,L,H)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     Subroutine E02BEV evaluates the (K+1) non-zero B-Splines of
C     degree K at T(L) <= X < T(L+1) using the stable recurrence
C     relation of De Boor and Cox.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           K, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  H(6), T(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F
      INTEGER           I, J
C     .. Local Arrays ..
      DOUBLE PRECISION  HH(5)
C     .. Executable Statements ..
C
      H(1) = ONE
      DO 60 J = 1, K
         DO 20 I = 1, J
            HH(I) = H(I)
   20    CONTINUE
         H(1) = ZERO
         DO 40 I = 1, J
            F = HH(I)/(T(L+I)-T(L+I-J))
            H(I) = H(I) + F*(T(L+I)-X)
            H(I+1) = F*(X-T(L+I-J))
   40    CONTINUE
   60 CONTINUE
      RETURN
      END
