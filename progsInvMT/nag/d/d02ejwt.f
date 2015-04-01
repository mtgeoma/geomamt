      DOUBLE PRECISION FUNCTION D02EJW(RDUM1,RDUM2)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 RDUM1, RDUM2
C     .. Scalars in Common ..
      CHARACTER*6                      CHDUM1, CHDUM2, GOPT
C     .. Common blocks ..
      COMMON                           /CD02EJ/CHDUM1, CHDUM2, GOPT
C     .. Save statement ..
      SAVE                             /CD02EJ/
C     .. Executable Statements ..
      D02EJW = 0.0D0
      GOPT = 'NOGOPT'
      RETURN
      END
