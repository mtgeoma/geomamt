      SUBROUTINE C06EAT(X0,PTS,X1,M1,X2,M2,X3,M3,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FOUR REAL FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C1, C2, C3, E, I1, I2, I3, IS0, IS1, IU0,
     *                  IU1, R1, R2, R3, RS0, RS1, RU0, RU1, S1, S2, S3,
     *                  TWOPI
      INTEGER           J, J0, J1, K, K0, K1, M4, MOVER2
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN, SQRT
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = (M-1)/2
      M4 = M*4
C
      DO 20 K = 1, PTS, M4
         RS0 = X0(K) + X1(K)
         RU0 = X0(K) - X1(K)
         RS1 = X2(K) + X3(K)
         RU1 = X2(K) - X3(K)
         X0(K) = RS0 + RS1
         X1(K) = RU0
         X2(K) = RS0 - RS1
         X3(K) = -RU1
   20 CONTINUE
C
      IF (MOVER2.LT.1) GO TO 80
      DO 60 J = 1, MOVER2
         J0 = J + 1
         J1 = M - 2*J
         ANGLE = TWOPI*DBLE(J)/DBLE(M4)
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
C
         DO 40 K0 = J0, PTS, M4
            K1 = K0 + J1
            R1 = X2(K0)*C1 + X2(K1)*S1
            I1 = X2(K1)*C1 - X2(K0)*S1
            R2 = X1(K0)*C2 + X1(K1)*S2
            I2 = X1(K1)*C2 - X1(K0)*S2
            R3 = X3(K0)*C3 + X3(K1)*S3
            I3 = X3(K1)*C3 - X3(K0)*S3
            RS0 = X0(K0) + R2
            IS0 = X0(K1) + I2
            RU0 = X0(K0) - R2
            IU0 = X0(K1) - I2
            RS1 = R1 + R3
            IS1 = I1 + I3
            RU1 = R1 - R3
            IU1 = I1 - I3
            X0(K0) = RS0 + RS1
            X3(K1) = IS0 + IS1
            X1(K0) = RU0 + IU1
            X2(K1) = IU0 - RU1
            X1(K1) = RS0 - RS1
            X2(K0) = -(IS0-IS1)
            X0(K1) = RU0 - IU1
            X3(K0) = -(IU0+RU1)
   40    CONTINUE
   60 CONTINUE
C
   80 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 120
C     E = COS(TWOPI/8.0)
      E = SQRT(0.5D0)
      J = M/2 + 1
      DO 100 K = J, PTS, M4
         RS0 = X0(K)
         RU0 = (X2(K)-X3(K))*E
         RS1 = (X2(K)+X3(K))*E
         RU1 = X1(K)
         X0(K) = RS0 + RU0
         X1(K) = RS0 - RU0
         X2(K) = RU1 - RS1
         X3(K) = -RU1 - RS1
  100 CONTINUE
C
  120 CONTINUE
      RETURN
      END
