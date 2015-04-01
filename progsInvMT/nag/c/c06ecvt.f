      SUBROUTINE C06ECV(X0,Y0,PTS,X1,Y1,M1,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX TWO COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), Y0(PTS), Y1(M1)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C, IS, IU, RS, RU, S, TWOPI
      INTEGER           J, K, K0, M2, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      M2 = M*2
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M2)
         ZERO = ANGLE .EQ. 0.0D0
         C = COS(ANGLE)
         S = SIN(ANGLE)
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         C = -C
   40    CONTINUE
C
         DO 100 K = K0, PTS, M2
            RS = X0(K) + X1(K)
            IS = Y0(K) + Y1(K)
            RU = X0(K) - X1(K)
            IU = Y0(K) - Y1(K)
            X0(K) = RS
            Y0(K) = IS
            IF (ZERO) GO TO 60
            X1(K) = RU*C + IU*S
            Y1(K) = IU*C - RU*S
            GO TO 80
   60       CONTINUE
            X1(K) = RU
            Y1(K) = IU
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
