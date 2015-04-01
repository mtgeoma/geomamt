      SUBROUTINE Y90PGF(LINE,VARCH)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==========================================
C         *  Y90PGF :  Print a Character Variable  *
C         ==========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      CHARACTER*(*)     LINE, VARCH
C     .. Local Scalars ..
      INTEGER           UNIT1
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
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
C
      IF (LEN(VARCH).LE.24) THEN
         REC(50:) = ' :  "'//VARCH(1:)//'"'
      ELSE
         REC(50:) = ' :  "'//VARCH(1:21)//'..."'
      END IF
C
      CALL X04BAF(UNIT1,REC)
C-----------------------------------------------------------------------
C
C     End of Y90PGF
C
C-----------------------------------------------------------------------
      RETURN
      END
