      SUBROUTINE Y90PAF(LINE)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ===================================
C         *  Y90PAF :  Print a Line of Text  *
C         ===================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
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
      CALL X04BAF(UNIT1,REC)
C-----------------------------------------------------------------------
C
C     End of Y90PAF
C
C-----------------------------------------------------------------------
      RETURN
      END
