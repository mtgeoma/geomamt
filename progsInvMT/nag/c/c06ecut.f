      SUBROUTINE C06ECU(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX THREE COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), Y0(PTS), Y1(M1), Y2(M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ANGLE, B, C1, C2, I0, I1, I2, IA, IB, IS, R0,
     *                  R1, R2, RA, RB, RS, S1, S2, T, TWOPI
      INTEGER           J, K, K0, M3, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN, SQRT
C     .. Executable Statements ..
      M3 = M*3
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
C     A = COS(TWOPI/3.0)
      A = -0.5D0
C     B = SIN(TWOPI/3.0)
      B = SQRT(0.75D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M3)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1*A + S1*B
         S1 = C1*B - S1*A
         C1 = T
         T = C2*A - S2*B
         S2 = -C2*B - S2*A
         C2 = T
   40    CONTINUE
C
         DO 100 K = K0, PTS, M3
            R0 = X0(K)
            I0 = Y0(K)
            RS = X1(K) + X2(K)
            IS = Y1(K) + Y2(K)
            X0(K) = R0 + RS
            Y0(K) = I0 + IS
            RA = R0 + RS*A
            IA = I0 + IS*A
            RB = (X1(K)-X2(K))*B
            IB = (Y1(K)-Y2(K))*B
            IF (ZERO) GO TO 60
            R1 = RA + IB
            I1 = IA - RB
            R2 = RA - IB
            I2 = IA + RB
            X1(K) = R1*C1 + I1*S1
            Y1(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            GO TO 80
   60       CONTINUE
            X1(K) = RA + IB
            Y1(K) = IA - RB
            X2(K) = RA - IB
            Y2(K) = IA + RB
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
