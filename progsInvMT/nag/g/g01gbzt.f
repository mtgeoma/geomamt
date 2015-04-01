      DOUBLE PRECISION FUNCTION G01GBZ(Z)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Calculates gamma(z+0.5)/gamma(z)
C     for z greater than 0.5
C
C     .. Parameters ..
      DOUBLE PRECISION                 HALF, BIG, SMALL
      PARAMETER                        (HALF=0.5D0,BIG=1.0D6,
     *                                 SMALL=25.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 Z
C     .. Local Scalars ..
      INTEGER                          IFAIL
C     .. External Functions ..
      DOUBLE PRECISION                 S14AAF
      EXTERNAL                         S14AAF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG
C     .. Executable Statements ..
C
C     for small z use call to S14AAF
C
      IF (Z.GE.0.5D0) THEN
         IF (Z.LE.SMALL) THEN
            IFAIL = 1
            G01GBZ = S14AAF(Z+HALF,IFAIL)/S14AAF(Z,IFAIL)
C
C     for large z use asymptotic formula
C
         ELSE IF (Z.LT.BIG) THEN
            G01GBZ = HALF*LOG(Z) + (1.0D0/(24.0D0*Z*Z)-1.0D0)/(8.0D0*Z)
            G01GBZ = EXP(G01GBZ)
         ELSE
            G01GBZ = HALF*LOG(Z) - 1.0D0/(8.0D0*Z)
            G01GBZ = EXP(G01GBZ)
         END IF
      ELSE
         G01GBZ = 0.0D0
      END IF
      END
