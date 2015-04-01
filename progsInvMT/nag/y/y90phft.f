      SUBROUTINE Y90PHF(LINE,IVAR)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ========================================
C         *  Y90PHF :  Print an Integer Variable  *
C         ========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVAR
      CHARACTER*(*)     LINE
C     .. Local Scalars ..
      INTEGER           UNIT1
      CHARACTER*80      REC
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
      WRITE (REC(54:),FMT='(I10)') IVAR
      CALL X04BAF(UNIT1,REC)
C-----------------------------------------------------------------------
C
C     End of Y90PHF
C
C-----------------------------------------------------------------------
      RETURN
      END
