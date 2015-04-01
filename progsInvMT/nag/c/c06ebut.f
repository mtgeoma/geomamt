      SUBROUTINE C06EBU(X0,PTS,X1,M1,X2,M2,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX THREE HERMITE FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ANGLE, B, C1, C2, I0, I1, I2, IA, IB, IS, R0,
     *                  R1, R2, RA, RB, RS, S1, S2, TWOPI
      INTEGER           J, J0, J1, K, K0, K1, M3, MOVER2
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN, SQRT
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = (M-1)/2
      M3 = M*3
C     A = COS(TWOPI/3.0)
      A = -0.5D0
C     B = 2.0*SIN(TWOPI/3.0)
      B = SQRT(3.0D0)
C
      DO 20 K = 1, PTS, M3
         R0 = X0(K)
         RS = 2.0D0*X1(K)
         RA = R0 + RS*A
         IB = X2(K)*B
         X0(K) = R0 + RS
         X1(K) = RA + IB
         X2(K) = RA - IB
   20 CONTINUE
C
      IF (MOVER2.LT.1) GO TO 80
      B = SIN(TWOPI/3.0D0)
      DO 60 J = 1, MOVER2
         J0 = J + 1
         J1 = M - 2*J
         ANGLE = TWOPI*DBLE(J)/DBLE(M3)
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
C
         DO 40 K0 = J0, PTS, M3
            K1 = K0 + J1
            R0 = X0(K0)
            I0 = X2(K1)
            RS = X1(K0) + X0(K1)
            IS = X1(K1) - X2(K0)
            RA = R0 + RS*A
            IA = I0 + IS*A
            RB = (X1(K0)-X0(K1))*B
            IB = (X1(K1)+X2(K0))*B
            X0(K0) = R0 + RS
            X0(K1) = I0 + IS
            R1 = RA + IB
            I1 = IA - RB
            R2 = RA - IB
            I2 = IA + RB
            X1(K0) = R1*C1 + I1*S1
            X1(K1) = I1*C1 - R1*S1
            X2(K0) = R2*C2 + I2*S2
            X2(K1) = I2*C2 - R2*S2
   40    CONTINUE
   60 CONTINUE
C
   80 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 120
C     A = COS(TWOPI/6.0)
      A = 0.5D0
C     B = 2.0*SIN(TWOPI/6.0)
      B = SQRT(3.0D0)
      J = M/2 + 1
      DO 100 K = J, PTS, M3
         R0 = X1(K)
         RS = 2.0D0*X0(K)
         RA = -R0 + RS*A
         IB = -X2(K)*B
         X0(K) = R0 + RS
         X1(K) = RA - IB
         X2(K) = -RA - IB
  100 CONTINUE
C
  120 CONTINUE
      RETURN
      END
