      SUBROUTINE Y90PKF(LINE,VARR)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ====================================
C         *  Y90PKF :  Print a Real Variable  *
C         ====================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  VARR
      CHARACTER*(*)     LINE
C     .. Local Scalars ..
      DOUBLE PRECISION  DUMMY
      INTEGER           UNIT1
      CHARACTER*80      REC
C     .. External Functions ..
      INTEGER           X02BEF
      EXTERNAL          X02BEF
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Print
C
C-----------------------------------------------------------------------
      CALL X04AAF(0,UNIT1)
C
      REC = ' '
      REC(5:) = LINE(1:)
      REC(50:) = ' :  '
      IF (X02BEF(DUMMY).LE.10) THEN
         WRITE (REC(54:),FMT='(1P,D13.6)') VARR
      ELSE
         WRITE (REC(54:),FMT='(1P,D19.12)') VARR
      END IF
      CALL X04BAF(UNIT1,REC)
C-----------------------------------------------------------------------
C
C     End of Y90PKF
C
C-----------------------------------------------------------------------
      RETURN
      END
