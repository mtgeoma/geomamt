      SUBROUTINE D02KDT(V,Y,PYP,K,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PYP, Y
      INTEGER           IFAIL, K
C     .. Array Arguments ..
      DOUBLE PRECISION  V(3)
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, PHI, R2
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, ATAN2, DBLE
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      B = V(1)
      IF (B.LE.ZER) GO TO 100
      R2 = Y*Y*B + PYP*PYP/B
      IF (R2.EQ.ZER) GO TO 80
      V(3) = LOG(R2)
      PHI = ATAN2(Y*B,PYP)
      IF (K.GE.0) GO TO 20
C     INITIAL BOUNDARY CONDITION
      IF (PHI.LT.ZER) PHI = PHI + PI
      IF (PHI.GE.PI) PHI = PHI - PI
      GO TO 40
C     FINAL BOUNDARY CONDITION
   20 IF (PHI.LE.ZER) PHI = PHI + PI
      IF (PHI.GT.PI) PHI = PHI - PI
      PHI = PHI + DBLE(K)*PI
   40 V(2) = TWO*PHI
      IFAIL = 0
   60 RETURN
   80 IFAIL = 1
      GO TO 60
  100 IFAIL = 2
      GO TO 60
      END
