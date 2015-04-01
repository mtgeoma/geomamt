      SUBROUTINE C06ECR(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,X4,Y4,M4,
     *                  X5,Y5,M5,X6,Y6,M6,X7,Y7,M7,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX EIGHT COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, M4, M5, M6, M7, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4), X5(M5),
     *                  X6(M6), X7(M7), Y0(PTS), Y1(M1), Y2(M2), Y3(M3),
     *                  Y4(M4), Y5(M5), Y6(M6), Y7(M7)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C1, C2, C3, C4, C5, C6, C7, E, I1, I2,
     *                  I3, I4, I5, I6, I7, IS0, IS1, IS2, IS3, ISS0,
     *                  ISS1, ISU0, ISU1, IU0, IU1, IU2, IU3, IUS0,
     *                  IUS1, IUU0, IUU1, R1, R2, R3, R4, R5, R6, R7,
     *                  RS0, RS1, RS2, RS3, RSS0, RSS1, RSU0, RSU1, RU0,
     *                  RU1, RU2, RU3, RUS0, RUS1, RUU0, RUU1, S1, S2,
     *                  S3, S4, S5, S6, S7, T, TWOPI
      INTEGER           J, K, K0, M8, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      M8 = M*8
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
      E = COS(TWOPI/8.0D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M8)
         ZERO = ANGLE .EQ. 0.0D0
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
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = (C1+S1)*E
         S1 = (C1-S1)*E
         C1 = T
         T = S2
         S2 = C2
         C2 = T
         T = (-C3+S3)*E
         S3 = (C3+S3)*E
         C3 = T
         C4 = -C4
         T = -(C5+S5)*E
         S5 = (-C5+S5)*E
         C5 = T
         T = -S6
         S6 = -C6
         C6 = T
         T = (C7-S7)*E
         S7 = -(C7+S7)*E
         C7 = T
   40    CONTINUE
C
         DO 100 K = K0, PTS, M8
            RS0 = X0(K) + X4(K)
            IS0 = Y0(K) + Y4(K)
            RU0 = X0(K) - X4(K)
            IU0 = Y0(K) - Y4(K)
            RS1 = X1(K) + X5(K)
            IS1 = Y1(K) + Y5(K)
            RU1 = X1(K) - X5(K)
            IU1 = Y1(K) - Y5(K)
            RS2 = X2(K) + X6(K)
            IS2 = Y2(K) + Y6(K)
            RU2 = X2(K) - X6(K)
            IU2 = Y2(K) - Y6(K)
            RS3 = X3(K) + X7(K)
            IS3 = Y3(K) + Y7(K)
            RU3 = X3(K) - X7(K)
            IU3 = Y3(K) - Y7(K)
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
            X0(K) = RSS0 + RSS1
            Y0(K) = ISS0 + ISS1
            IF (ZERO) GO TO 60
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
            X4(K) = R1*C1 + I1*S1
            Y4(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            X6(K) = R3*C3 + I3*S3
            Y6(K) = I3*C3 - R3*S3
            X1(K) = R4*C4 + I4*S4
            Y1(K) = I4*C4 - R4*S4
            X5(K) = R5*C5 + I5*S5
            Y5(K) = I5*C5 - R5*S5
            X3(K) = R6*C6 + I6*S6
            Y3(K) = I6*C6 - R6*S6
            X7(K) = R7*C7 + I7*S7
            Y7(K) = I7*C7 - R7*S7
            GO TO 80
   60       CONTINUE
            X4(K) = RUU0 + RUU1
            Y4(K) = IUU0 + IUU1
            X2(K) = RSU0 + ISU1
            Y2(K) = ISU0 - RSU1
            X6(K) = RUS0 + IUS1
            Y6(K) = IUS0 - RUS1
            X1(K) = RSS0 - RSS1
            Y1(K) = ISS0 - ISS1
            X5(K) = RUU0 - RUU1
            Y5(K) = IUU0 - IUU1
            X3(K) = RSU0 - ISU1
            Y3(K) = ISU0 + RSU1
            X7(K) = RUS0 - IUS1
            Y7(K) = IUS0 + RUS1
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
