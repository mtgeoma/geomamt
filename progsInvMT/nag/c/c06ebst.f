      SUBROUTINE C06EBS(X0,PTS,X1,M1,X2,M2,X3,M3,X4,M4,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FIVE HERMITE FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, M4, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, ANGLE, AS, AU, B1, B2, C1, C2, C3, C4,
     *                  I0, I1, I2, I3, I4, IA1, IA2, IAS, IAU, IB1,
     *                  IB2, IS1, IS2, ISS, IU1, IU2, R0, R1, R2, R3,
     *                  R4, RA1, RA2, RAS, RAU, RB1, RB2, RS1, RS2, RSS,
     *                  RU1, RU2, S1, S2, S3, S4, TWOPI
      INTEGER           J, J0, J1, K, K0, K1, M5, MOVER2
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN, SQRT
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = (M-1)/2
      M5 = M*5
      A1 = COS(TWOPI/5.0D0)
      B1 = SIN(TWOPI/5.0D0)
      A2 = COS(2.0D0*TWOPI/5.0D0)
      B2 = SIN(2.0D0*TWOPI/5.0D0)
      AS = -1.0D0/4.0D0
      AU = SQRT(5.0D0)/4.0D0
C
      DO 20 K = 1, PTS, M5
         R0 = X0(K)
         RS1 = 2.0D0*X1(K)
         IU1 = 2.0D0*X4(K)
         RS2 = 2.0D0*X2(K)
         IU2 = 2.0D0*X3(K)
         RSS = RS1 + RS2
         RAS = R0 + RSS*AS
         RAU = (RS1-RS2)*AU
         RA1 = RAS + RAU
         RA2 = RAS - RAU
         IB1 = IU1*B1 + IU2*B2
         IB2 = IU1*B2 - IU2*B1
         X0(K) = R0 + RSS
         X1(K) = RA1 + IB1
         X4(K) = RA1 - IB1
         X2(K) = RA2 + IB2
         X3(K) = RA2 - IB2
   20 CONTINUE
C
      IF (MOVER2.LT.1) GO TO 80
      DO 60 J = 1, MOVER2
         J0 = J + 1
         J1 = M - 2*J
         ANGLE = TWOPI*DBLE(J)/DBLE(M5)
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         C4 = C2*C2 - S2*S2
         S4 = S2*C2 + C2*S2
C
         DO 40 K0 = J0, PTS, M5
            K1 = K0 + J1
            R0 = X0(K0)
            I0 = X4(K1)
            RS1 = X1(K0) + X0(K1)
            IS1 = X3(K1) - X4(K0)
            RU1 = X1(K0) - X0(K1)
            IU1 = X3(K1) + X4(K0)
            RS2 = X2(K0) + X1(K1)
            IS2 = X2(K1) - X3(K0)
            RU2 = X2(K0) - X1(K1)
            IU2 = X2(K1) + X3(K0)
            RSS = RS1 + RS2
            ISS = IS1 + IS2
            RAS = R0 + RSS*AS
            IAS = I0 + ISS*AS
            RAU = (RS1-RS2)*AU
            IAU = (IS1-IS2)*AU
            RA1 = RAS + RAU
            IA1 = IAS + IAU
            RA2 = RAS - RAU
            IA2 = IAS - IAU
            RB1 = RU1*B1 + RU2*B2
            IB1 = IU1*B1 + IU2*B2
            RB2 = RU1*B2 - RU2*B1
            IB2 = IU1*B2 - IU2*B1
            X0(K0) = R0 + RSS
            X0(K1) = I0 + ISS
            R1 = RA1 + IB1
            I1 = IA1 - RB1
            R2 = RA2 + IB2
            I2 = IA2 - RB2
            R3 = RA2 - IB2
            I3 = IA2 + RB2
            R4 = RA1 - IB1
            I4 = IA1 + RB1
            X1(K0) = R1*C1 + I1*S1
            X1(K1) = I1*C1 - R1*S1
            X2(K0) = R2*C2 + I2*S2
            X2(K1) = I2*C2 - R2*S2
            X3(K0) = R3*C3 + I3*S3
            X3(K1) = I3*C3 - R3*S3
            X4(K0) = R4*C4 + I4*S4
            X4(K1) = I4*C4 - R4*S4
   40    CONTINUE
   60 CONTINUE
C
   80 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 120
      A1 = COS(TWOPI/10.0D0)
      B1 = SIN(TWOPI/10.0D0)
      A2 = COS(TWOPI/5.0D0)
      B2 = SIN(TWOPI/5.0D0)
      J = M/2 + 1
      DO 100 K = J, PTS, M5
         R0 = X2(K)
         RS1 = 2.0D0*X0(K)
         IU1 = -2.0D0*X4(K)
         RS2 = 2.0D0*X1(K)
         IU2 = -2.0D0*X3(K)
         RA1 = -R0 + RS1*A1 - RS2*A2
         RA2 = R0 + RS1*A2 - RS2*A1
         IB1 = IU1*B1 + IU2*B2
         IB2 = IU1*B2 - IU2*B1
         X0(K) = R0 + RS1 + RS2
         X1(K) = RA1 - IB1
         X4(K) = -RA1 - IB1
         X2(K) = RA2 - IB2
         X3(K) = -RA2 - IB2
  100 CONTINUE
C
  120 CONTINUE
      RETURN
      END
