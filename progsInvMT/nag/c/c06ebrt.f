      SUBROUTINE C06EBR(X0,PTS,X1,M1,X2,M2,X3,M3,X4,M4,X5,M5,X6,M6,X7,
     *                  M7,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX EIGHT HERMITE FOURIER TRANSFORM KERNEL
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
         RS0 = X0(K) + X4(K)
         RU0 = X0(K) - X4(K)
         RS1 = X1(K) + X3(K)
         IS1 = X7(K) - X5(K)
         RU1 = X1(K) - X3(K)
         IU1 = X7(K) + X5(K)
         RS2 = 2.0D0*X2(K)
         IU2 = 2.0D0*X6(K)
         RSS0 = RS0 + RS2
         RSU0 = RS0 - RS2
         RSS1 = 2.0D0*RS1
         ISU1 = 2.0D0*IS1
         RUS0 = RU0 - IU2
         RUU0 = RU0 + IU2
         IUS1 = -2.0D0*E*(RU1-IU1)
         RUU1 = 2.0D0*E*(RU1+IU1)
         X0(K) = RSS0 + RSS1
         X4(K) = RUU0 + RUU1
         X2(K) = RSU0 + ISU1
         X6(K) = RUS0 + IUS1
         X1(K) = RSS0 - RSS1
         X5(K) = RUU0 - RUU1
         X3(K) = RSU0 - ISU1
         X7(K) = RUS0 - IUS1
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
            RS0 = X0(K0) + X3(K1)
            IS0 = X7(K1) - X4(K0)
            RU0 = X0(K0) - X3(K1)
            IU0 = X7(K1) + X4(K0)
            RS1 = X1(K0) + X2(K1)
            IS1 = X6(K1) - X5(K0)
            RU1 = X1(K0) - X2(K1)
            IU1 = X6(K1) + X5(K0)
            RS2 = X2(K0) + X1(K1)
            IS2 = X5(K1) - X6(K0)
            RU2 = X2(K0) - X1(K1)
            IU2 = X5(K1) + X6(K0)
            RS3 = X3(K0) + X0(K1)
            IS3 = X4(K1) - X7(K0)
            RU3 = X3(K0) - X0(K1)
            IU3 = X4(K1) + X7(K0)
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
            X0(K1) = ISS0 + ISS1
            R1 = RUU0 + RUU1
            I1 = IUU0 + IUU1
            R2 = RSU0 + ISU1
            I2 = ISU0 - RSU1
            R3 = RUS0 + IUS1
            I3 = IUS0 - RUS1
            R4 = RSS0 - RSS1
            I4 = ISS0 - ISS1
            R5 = RUU0 - RUU1
            I5 = IUU0 - IUU1
            R6 = RSU0 - ISU1
            I6 = ISU0 + RSU1
            R7 = RUS0 - IUS1
            I7 = IUS0 + RUS1
            X4(K0) = R1*C1 + I1*S1
            X4(K1) = I1*C1 - R1*S1
            X2(K0) = R2*C2 + I2*S2
            X2(K1) = I2*C2 - R2*S2
            X6(K0) = R3*C3 + I3*S3
            X6(K1) = I3*C3 - R3*S3
            X1(K0) = R4*C4 + I4*S4
            X1(K1) = I4*C4 - R4*S4
            X5(K0) = R5*C5 + I5*S5
            X5(K1) = I5*C5 - R5*S5
            X3(K0) = R6*C6 + I6*S6
            X3(K1) = I6*C6 - R6*S6
            X7(K0) = R7*C7 + I7*S7
            X7(K1) = I7*C7 - R7*S7
   40    CONTINUE
   60 CONTINUE
C
   80 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 120
      A = COS(TWOPI/16.0D0)
      B = COS(3.0D0*TWOPI/16.0D0)
      J = M/2 + 1
      DO 100 K = J, PTS, M8
         RS0 = 2.0D0*(X0(K)+X3(K))
         IS0 = 2.0D0*(X7(K)+X4(K))
         RU0 = 2.0D0*(X0(K)-X3(K))
         IU0 = 2.0D0*(X7(K)-X4(K))
         RS1 = 2.0D0*(X1(K)+X2(K))
         IS1 = 2.0D0*(X6(K)+X5(K))
         RU1 = 2.0D0*(X1(K)-X2(K))
         IU1 = 2.0D0*(X6(K)-X5(K))
         X0(K) = RS0 + RS1
         X1(K) = IU0 - IU1
         ISS0 = IS0*B + IS1*A
         ISU0 = IS0*A - IS1*B
         RSU1 = (RS0-RS1)*E
         RUS0 = RU0*A + RU1*B
         RUU0 = RU0*B - RU1*A
         IUS1 = (IU0+IU1)*E
         X4(K) = RUS0 + ISS0
         X2(K) = RSU1 + IUS1
         X6(K) = RUU0 + ISU0
         X5(K) = -RUU0 + ISU0
         X3(K) = -RSU1 + IUS1
         X7(K) = -RUS0 + ISS0
  100 CONTINUE
C
  120 CONTINUE
      RETURN
      END
