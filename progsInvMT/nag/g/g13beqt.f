      SUBROUTINE G13BEQ(F,ALPHA,R,NPD,V)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEQ CALCULATES VALUES OF RHO AND V WHICH
C     ARE USED IN AUTOREGRESSIVE PARAMETER DERIVATIVE CORRECTIONS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  V
      INTEGER           NPD
C     .. Array Arguments ..
      DOUBLE PRECISION  ALPHA(NPD), F(NPD), R(NPD)
C     .. Local Scalars ..
      DOUBLE PRECISION  TAU, U
      INTEGER           I, J, K, KM
C     .. Data statements ..
      DATA              U/1.0D0/
C     .. Executable Statements ..
      V = U
      DO 80 K = 1, NPD
         TAU = (U-F(K))*(U+F(K))
         R(K) = (-V)*F(K)
         V = V*TAU
         IF (K.EQ.1) GO TO 80
         KM = K - 1
         DO 20 J = 1, KM
            I = K - J
            R(K) = R(K) - F(J)*R(I)
   20    CONTINUE
         DO 40 J = 1, KM
            I = K - J
            ALPHA(J) = F(J) + F(K)*F(I)
   40    CONTINUE
         DO 60 J = 1, KM
            F(J) = ALPHA(J)
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
