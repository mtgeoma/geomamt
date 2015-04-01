      SUBROUTINE C06EAR(X0,PTS,X1,M1,X2,M2,X3,M3,X4,M4,X5,M5,X6,M6,X7,
     *                  M7,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX EIGHT REAL FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, M4, M5, M6, M7, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4), X5(M5),
     *                  X6(M6), X7(M7)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ANGLE, B, C1, C2, C3, C4, C5, C6, C7, E, I1,
     *                  I2, I3, I4, I5, I6, I7, IS0, IS1, IS2, IS3,
     *                  ISS0, ISS1, ISU0, ISU1, IU0, IU1, IU2, IU3,
     *                  IUS0, IUS1, IUU0, IUU1, R1, R2, R3, R4, R5, R6,
     *                  R7, RS0, RS1, RS2, RS3, RSS0, RSS1, RSU0, RSU1,
     *                  RU0, RU1, RU2, RU3, RUS0, RUS1, RUU0, RUU1, S1,
     *                  S2, S3, S4, S5, S6, S7, T, TWOPI
      INTEGER           J, J0, J1, K, K0, K1, M8, MOVER2
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN, SQRT
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = (M-1)/2
      M8 = M*8
C     E = COS(TWOPI/8.0)
      E = SQRT(0.5D0)
C
      DO 20 K = 1, PTS, M8
         RS0 = X0(K) + X1(K)
         RU0 = X0(K) - X1(K)
         RS1 = X4(K) + X5(K)
         RU1 = X4(K) - X5(K)
         RS2 = X2(K) + X3(K)
         RU2 = X2(K) - X3(K)
         RS3 = X6(K) + X7(K)
         RU3 = X6(K) - X7(K)
         RSS0 = RS0 + RS2
         RSU0 = RS0 - RS2
         RSS1 = RS1 + RS3
         RSU1 = RS1 - RS3
         RUU0 = RU0
         IUU0 = -RU2
         RUU1 = (RU1-RU3)*E
         IUU1 = -(RU1+RU3)*E
         X0(K) = RSS0 + RSS1
         X1(K) = RUU0 + RUU1
         X2(K) = RSU0
         X3(K) = RUU0 - RUU1
         X4(K) = RSS0 - RSS1
         X5(K) = -IUU0 + IUU1
         X6(K) = -RSU1
         X7(K) = IUU0 + IUU1
   20 CONTINUE
C
      IF (MOVER2.LT.1) GO TO 80
      DO 60 J = 1, MOVER2
         J0 = J + 1
         J1 = M - 2*J
         ANGLE = TWOPI*DBLE(J)/DBLE(M8)
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         C4 = C2*C2 - S2*S2
         S4 = S2*C2 + C2*S2
         C5 = C4*C1 - S4*S1
         S5 = S4*C1 + C4*S1
         C6 = C4*C2 - S4*S2
         S6 = S4*C2 + C4*S2
         C7 = C4*C3 - S4*S3
         S7 = S4*C3 + C4*S3
C
         DO 40 K0 = J0, PTS, M8
            K1 = K0 + J1
            R1 = X4(K0)*C1 + X4(K1)*S1
            I1 = X4(K1)*C1 - X4(K0)*S1
            R2 = X2(K0)*C2 + X2(K1)*S2
            I2 = X2(K1)*C2 - X2(K0)*S2
            R3 = X6(K0)*C3 + X6(K1)*S3
            I3 = X6(K1)*C3 - X6(K0)*S3
            R4 = X1(K0)*C4 + X1(K1)*S4
            I4 = X1(K1)*C4 - X1(K0)*S4
            R5 = X5(K0)*C5 + X5(K1)*S5
            I5 = X5(K1)*C5 - X5(K0)*S5
            R6 = X3(K0)*C6 + X3(K1)*S6
            I6 = X3(K1)*C6 - X3(K0)*S6
            R7 = X7(K0)*C7 + X7(K1)*S7
            I7 = X7(K1)*C7 - X7(K0)*S7
            RS0 = X0(K0) + R4
            IS0 = X0(K1) + I4
            RU0 = X0(K0) - R4
            IU0 = X0(K1) - I4
            RS1 = R1 + R5
            IS1 = I1 + I5
            RU1 = R1 - R5
            IU1 = I1 - I5
            RS2 = R2 + R6
            IS2 = I2 + I6
            RU2 = R2 - R6
            IU2 = I2 - I6
            RS3 = R3 + R7
            IS3 = I3 + I7
            RU3 = R3 - R7
            IU3 = I3 - I7
            RSS0 = RS0 + RS2
            ISS0 = IS0 + IS2
            RSU0 = RS0 - RS2
            ISU0 = IS0 - IS2
            RSS1 = RS1 + RS3
            ISS1 = IS1 + IS3
            RSU1 = RS1 - RS3
            ISU1 = IS1 - IS3
            RUS0 = RU0 - IU2
            IUS0 = IU0 + RU2
            RUU0 = RU0 + IU2
            IUU0 = IU0 - RU2
            RUS1 = RU1 - IU3
            IUS1 = IU1 + RU3
            RUU1 = RU1 + IU3
            IUU1 = IU1 - RU3
            T = (RUS1+IUS1)*E
            IUS1 = (IUS1-RUS1)*E
            RUS1 = T
            T = (RUU1+IUU1)*E
            IUU1 = (IUU1-RUU1)*E
            RUU1 = T
            X0(K0) = RSS0 + RSS1
            X7(K1) = ISS0 + ISS1
            X1(K0) = RUU0 + RUU1
            X6(K1) = IUU0 + IUU1
            X2(K0) = RSU0 + ISU1
            X5(K1) = ISU0 - RSU1
            X3(K0) = RUS0 + IUS1
            X4(K1) = IUS0 - RUS1
            X3(K1) = RSS0 - RSS1
            X4(K0) = -(ISS0-ISS1)
            X2(K1) = RUU0 - RUU1
            X5(K0) = -(IUU0-IUU1)
            X1(K1) = RSU0 - ISU1
            X6(K0) = -(ISU0+RSU1)
            X0(K1) = RUS0 - IUS1
            X7(K0) = -(IUS0+RUS1)
   40    CONTINUE
   60 CONTINUE
C
   80 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 120
      A = COS(TWOPI/16.0D0)
      B = COS(3.0D0*TWOPI/16.0D0)
      J = M/2 + 1
      DO 100 K = J, PTS, M8
         RS1 = X4(K) + X7(K)
         RU1 = X4(K) - X7(K)
         RS2 = (X2(K)+X3(K))*E
         RU2 = (X2(K)-X3(K))*E
         RS3 = X6(K) + X5(K)
         RU3 = X6(K) - X5(K)
         RSS0 = X1(K) + RS2
         RSU0 = X1(K) - RS2
         RSS1 = RS1*B + RS3*A
         RSU1 = RS1*A - RS3*B
         RUS0 = X0(K) + RU2
         RUU0 = X0(K) - RU2
         RUS1 = RU1*A + RU3*B
         RUU1 = RU1*B - RU3*A
         X0(K) = RUS0 + RUS1
         X1(K) = RUU0 + RUU1
         X2(K) = RUU0 - RUU1
         X3(K) = RUS0 - RUS1
         X4(K) = RSS0 - RSS1
         X5(K) = -RSU0 - RSU1
         X6(K) = RSU0 - RSU1
         X7(K) = -RSS0 - RSS1
  100 CONTINUE
C
  120 CONTINUE
      RETURN
      END
