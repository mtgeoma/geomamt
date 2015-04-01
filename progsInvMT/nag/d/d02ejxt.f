      SUBROUTINE D02EJX(RDUM1,RDUM2)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RDUM1, RDUM2
C     .. Scalars in Common ..
      CHARACTER*6       CHDUM1, CHDUM2, OUTOPT
C     .. Common blocks ..
      COMMON            /CD02EJ/CHDUM1, OUTOPT, CHDUM2
C     .. Save statement ..
      SAVE              /CD02EJ/
C     .. Executable Statements ..
      OUTOPT = 'NOXOUT'
      RETURN
      END
