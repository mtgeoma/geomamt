      SUBROUTINE E02GBV(L,M,N,E,IER,MPL1,RES,GRD,PEN,EL1N,PENPAR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     CL1 VERSION OF E02GBV.
C
C     THIS ROUTINE ADMINISTERS THE EVALUATION OF THE
C     PENALTY (OBJECTIVE) FUNCTION GIVEN THE EQUATION
C     AND CONSTRAINT RESIDUALS.  IT ALSO COMPUTES THE
C     RESTRICTED GRADIENT OF THE FUNCTION.
C     ***************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EL1N, PEN, PENPAR
      INTEGER           IER, L, M, MPL1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  E(IER,MPL1), GRD(N), RES(MPL1)
C     .. Local Scalars ..
      DOUBLE PRECISION  DKSTMP, ONE, ZERO
      INTEGER           I, MP1, MPL
C     .. External Subroutines ..
      EXTERNAL          E02GBL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN
C     .. Data statements ..
      DATA              ZERO/0.0D+00/, ONE/1.0D+00/
C     .. Executable Statements ..
      MPL = M + L
      MP1 = M + 1
      EL1N = ZERO
      DO 20 I = 1, N
         GRD(I) = ZERO
   20 CONTINUE
      IF (1.GT.M) GO TO 60
      DO 40 I = 1, M
         IF (RES(I).EQ.ZERO) GO TO 40
         EL1N = EL1N + ABS(RES(I))
         DKSTMP = 1.0D0
         IF (RES(I).NE.0.0D0) DKSTMP = SIGN(ONE,RES(I))
         CALL E02GBL(N,DKSTMP*PENPAR,E(1,I),GRD)
   40 CONTINUE
   60 CONTINUE
      PEN = PENPAR*EL1N
      IF (MP1.GT.MPL) GO TO 100
      DO 80 I = MP1, MPL
         IF (RES(I).GE.ZERO) GO TO 80
         PEN = PEN + ABS(RES(I))
         CALL E02GBL(N,-ONE,E(1,I),GRD)
   80 CONTINUE
  100 CONTINUE
      RETURN
      END
