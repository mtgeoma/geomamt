      SUBROUTINE D02EJY(RDUM1,RDUM2,RDUM3)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RDUM1
C     .. Array Arguments ..
      DOUBLE PRECISION  RDUM2(*), RDUM3(*)
C     .. Scalars in Common ..
      CHARACTER*6       CHDUM1, CHDUM2, JACOPT
C     .. Common blocks ..
      COMMON            /CD02EJ/JACOPT, CHDUM1, CHDUM2
C     .. Save statement ..
      SAVE              /CD02EJ/
C     .. Executable Statements ..
      JACOPT = 'NUMERI'
      RETURN
      END
