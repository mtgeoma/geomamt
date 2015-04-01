      SUBROUTINE C06EBT(X0,PTS,X1,M1,X2,M2,X3,M3,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FOUR HERMITE FOURIER TRANSFORM KERNEL
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
         RS0 = X0(K) + X2(K)
         RU0 = X0(K) - X2(K)
         RS1 = 2.0D0*X1(K)
         IU1 = 2.0D0*X3(K)
         X0(K) = RS0 + RS1
         X2(K) = RU0 + IU1
         X1(K) = RS0 - RS1
         X3(K) = RU0 - IU1
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
            RS0 = X0(K0) + X1(K1)
            IS0 = X3(K1) - X2(K0)
            RU0 = X0(K0) - X1(K1)
            IU0 = X3(K1) + X2(K0)
            RS1 = X1(K0) + X0(K1)
            IS1 = X2(K1) - X3(K0)
            RU1 = X1(K0) - X0(K1)
            IU1 = X2(K1) + X3(K0)
            X0(K0) = RS0 + RS1
            X0(K1) = IS0 + IS1
            R1 = RU0 + IU1
            I1 = IU0 - RU1
            R2 = RS0 - RS1
            I2 = IS0 - IS1
            R3 = RU0 - IU1
            I3 = IU0 + RU1
            X2(K0) = R1*C1 + I1*S1
            X2(K1) = I1*C1 - R1*S1
            X1(K0) = R2*C2 + I2*S2
            X1(K1) = I2*C2 - R2*S2
            X3(K0) = R3*C3 + I3*S3
            X3(K1) = I3*C3 - R3*S3
   40    CONTINUE
   60 CONTINUE
C
   80 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 120
C     E = 2.0*COS(TWOPI/8.0)
      E = SQRT(2.0D0)
      J = M/2 + 1
      DO 100 K = J, PTS, M4
         RS0 = 2.0D0*(X0(K)+X1(K))
         IS0 = (X3(K)+X2(K))*E
         RU0 = (X0(K)-X1(K))*E
         IU0 = 2.0D0*(X3(K)-X2(K))
         X0(K) = RS0
         X2(K) = RU0 + IS0
         X1(K) = IU0
         X3(K) = -RU0 + IS0
  100 CONTINUE
C
  120 CONTINUE
      RETURN
      END
