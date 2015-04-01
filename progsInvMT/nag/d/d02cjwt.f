      DOUBLE PRECISION FUNCTION D02CJW(RDUM1,RDUM2)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO
      PARAMETER                        (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 RDUM1, RDUM2
C     .. Scalars in Common ..
      CHARACTER*6                      CHDUM, GOPT
C     .. Common blocks ..
      COMMON                           /AD02CJ/CHDUM, GOPT
C     .. Save statement ..
      SAVE                             /AD02CJ/
C     .. Executable Statements ..
      D02CJW = ZERO
      GOPT = 'NOGOPT'
      RETURN
      END
