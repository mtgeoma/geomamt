      SUBROUTINE D02KDV(Y,BNEW,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-499 (AUG 1986).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BNEW
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(3)
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, BUP, CPHI, SPHI
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, ATAN2, COS, SIN
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      B = Y(1)
      IF (B.LE.ZER .OR. BNEW.LE.ZER) GO TO 20
      Y(1) = BNEW
      BUP = B/BNEW
      B = BNEW/B
      CPHI = COS(Y(2))
      SPHI = SIN(Y(2))
      Y(2) = Y(2) + TWO*ATAN2((B-ONE)*SPHI,B+ONE-(B-ONE)*CPHI)
      Y(3) = Y(3) + LOG((B+BUP-(B-BUP)*CPHI)/TWO)
      IFAIL = 0
      RETURN
   20 IFAIL = 1
      RETURN
      END
