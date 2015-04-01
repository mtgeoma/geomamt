      SUBROUTINE E02BEW(T,N,K2,B,NEST)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     Subroutine E02BEW calculates the discontinuity jumps of the Kth
C     derivative of the B-Splines of degree K at the knots
C     T(K+2)..T(N-K-1)
C     .. Scalar Arguments ..
      INTEGER           K2, N, NEST
C     .. Array Arguments ..
      DOUBLE PRECISION  B(NEST,K2), T(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  FAC, PROD
      INTEGER           I, J, K, K1, L, NRINT
C     .. Local Arrays ..
      DOUBLE PRECISION  H(12)
C     .. Executable Statements ..
C
      K1 = K2 - 1
      K = K1 - 1
      NRINT = N - K1 - K
      FAC = NRINT/(T(N-K1+1)-T(K1))
      DO 80 L = K2, N - K1
         DO 20 J = 1, K1
            H(J) = T(L) - T(L+J-K2)
            H(J+K1) = T(L) - T(L+J)
   20    CONTINUE
         DO 60 J = 1, K2
            PROD = H(J)
            DO 40 I = J + 1, J + K
               PROD = PROD*H(I)*FAC
   40       CONTINUE
            B(L-K1,J) = (T(L+J-1)-T(L-K1+J-1))/PROD
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
