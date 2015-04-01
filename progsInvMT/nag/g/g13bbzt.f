      SUBROUTINE G13BBZ(MR,PAR,NPAR,EF,IERR)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C        ROUTINE TO CHECK TF MODEL VALIDITY
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EF
      INTEGER           IERR, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  PAR(NPAR)
      INTEGER           MR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, PG, Q
      INTEGER           K, KC, L
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          G13AEX
C     .. Executable Statements ..
      IERR = 0
      Q = 1.0D0
      EPS = EF*X02AJF()
      K = MR(3)
      L = MR(2) + 2
      CALL G13AEX(PAR(L),K,EPS,PG,KC)
      IF (KC.GT.0) GO TO 20
      IERR = -1
   20 RETURN
      END
