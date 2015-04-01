      SUBROUTINE G13BAY(MR,PAR,NPAR,EF,IERR)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     G13BAY CHECKS ARIMA MODEL VALIDITY FOR G13BAF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EF
      INTEGER           IERR, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  PAR(NPAR)
      INTEGER           MR(7)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, PG
      INTEGER           I, K, KC, L1
C     .. Local Arrays ..
      INTEGER           N(4)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          G13AEX
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = 0
      EPS = EF*X02AJF()
      N(1) = 1
      N(2) = 3
      N(3) = 4
      N(4) = 6
      L1 = 1
      DO 20 I = 1, 4
         K = N(I)
         K = MR(K)
         IF (K.EQ.0) GO TO 20
         CALL G13AEX(PAR(L1),K,EPS,PG,KC)
         IF (MOD(I,2).EQ.0 .AND. KC.EQ.0 .OR. KC.LT.0) GO TO 40
         L1 = L1 + K
   20 CONTINUE
      GO TO 60
   40 IERR = -1
   60 RETURN
      END
