      SUBROUTINE D02CJX(RDUM1,RDUM2)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RDUM1, RDUM2
C     .. Scalars in Common ..
      CHARACTER*6       CHDUM, OUTOPT
C     .. Common blocks ..
      COMMON            /AD02CJ/OUTOPT, CHDUM
C     .. Save statement ..
      SAVE              /AD02CJ/
C     .. Executable Statements ..
      OUTOPT = 'NOXOUT'
      RETURN
      END
