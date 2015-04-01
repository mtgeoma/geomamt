      SUBROUTINE C06EBV(X0,PTS,X1,M1,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX TWO HERMITE FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C, IS, IU, RS, RU, S, TWOPI
      INTEGER           J, J0, J1, K, K0, K1, M2, MOVER2
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = (M-1)/2
      M2 = M*2
C
      DO 20 K = 1, PTS, M2
         RS = X0(K) + X1(K)
         RU = X0(K) - X1(K)
         X0(K) = RS
         X1(K) = RU
   20 CONTINUE
C
      IF (MOVER2.LT.1) GO TO 80
      DO 60 J = 1, MOVER2
         J0 = J + 1
         J1 = M - 2*J
         ANGLE = TWOPI*DBLE(J)/DBLE(M2)
         C = COS(ANGLE)
         S = SIN(ANGLE)
C
         DO 40 K0 = J0, PTS, M2
            K1 = K0 + J1
            RS = X0(K0) + X0(K1)
            IS = X1(K1) - X1(K0)
            RU = X0(K0) - X0(K1)
            IU = X1(K1) + X1(K0)
            X0(K0) = RS
            X0(K1) = IS
            X1(K0) = RU*C + IU*S
            X1(K1) = IU*C - RU*S
   40    CONTINUE
   60 CONTINUE
C
   80 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 120
      J = M/2 + 1
      DO 100 K = J, PTS, M2
         X0(K) = 2.0D0*X0(K)
         X1(K) = 2.0D0*X1(K)
  100 CONTINUE
C
  120 CONTINUE
      RETURN
      END
