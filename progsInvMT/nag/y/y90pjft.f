      SUBROUTINE Y90PJF(LINE,LVAR)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90PJF :  Print a Logical Variable  *
C         =======================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      LOGICAL           LVAR
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
      IF (LVAR) THEN
         REC(54:) = ' .TRUE.'
      ELSE
         REC(54:) = ' .FALSE.'
      END IF
      CALL X04BAF(UNIT1,REC)
C-----------------------------------------------------------------------
C
C     End of Y90PJF
C
C-----------------------------------------------------------------------
      RETURN
      END
