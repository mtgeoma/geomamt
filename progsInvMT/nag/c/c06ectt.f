      SUBROUTINE C06ECT(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FOUR COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), Y0(PTS),
     *                  Y1(M1), Y2(M2), Y3(M3)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C1, C2, C3, I1, I2, I3, IS0, IS1, IU0,
     *                  IU1, R1, R2, R3, RS0, RS1, RU0, RU1, S1, S2, S3,
     *                  T, TWOPI
      INTEGER           J, K, K0, M4, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      M4 = M*4
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M4)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1
         C1 = S1
         S1 = T
         C2 = -C2
         T = C3
         C3 = -S3
         S3 = -T
   40    CONTINUE
C
         DO 100 K = K0, PTS, M4
            RS0 = X0(K) + X2(K)
            IS0 = Y0(K) + Y2(K)
            RU0 = X0(K) - X2(K)
            IU0 = Y0(K) - Y2(K)
            RS1 = X1(K) + X3(K)
            IS1 = Y1(K) + Y3(K)
            RU1 = X1(K) - X3(K)
            IU1 = Y1(K) - Y3(K)
            X0(K) = RS0 + RS1
            Y0(K) = IS0 + IS1
            IF (ZERO) GO TO 60
            R1 = RU0 + IU1
            I1 = IU0 - RU1
            R2 = RS0 - RS1
            I2 = IS0 - IS1
            R3 = RU0 - IU1
            I3 = IU0 + RU1
            X2(K) = R1*C1 + I1*S1
            Y2(K) = I1*C1 - R1*S1
            X1(K) = R2*C2 + I2*S2
            Y1(K) = I2*C2 - R2*S2
            X3(K) = R3*C3 + I3*S3
            Y3(K) = I3*C3 - R3*S3
            GO TO 80
   60       CONTINUE
            X2(K) = RU0 + IU1
            Y2(K) = IU0 - RU1
            X1(K) = RS0 - RS1
            Y1(K) = IS0 - IS1
            X3(K) = RU0 - IU1
            Y3(K) = IU0 + RU1
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
