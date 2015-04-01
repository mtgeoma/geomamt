      SUBROUTINE Y90HAF(INAG,NAG,TNAG,WNAG,XNAG)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =============================================
C         *  Y90HAF :  Print message for NAG testing  *
C         =============================================
C
C
C     Print a NAG routine related message
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           INAG, XNAG
      LOGICAL           TNAG, WNAG
      CHARACTER*6       NAG
C     .. Local Scalars ..
      INTEGER           UNIT1
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Print message
C
C-----------------------------------------------------------------------
      CALL X04AAF(0,UNIT1)
C
      IF (XNAG.LE.1) THEN
         IF (TNAG) THEN
            WRITE (REC,FMT=99999) INAG, NAG
         ELSE
            WRITE (REC,FMT=99998) INAG, NAG
         END IF
         CALL X04BAF(UNIT1,REC)
         IF ( .NOT. WNAG) THEN
            REC = ' '
            REC(19:) = '++  WARNINGS HAVE BEEN RAISED'
            CALL X04BAF(UNIT1,REC)
         END IF
      ELSE IF (XNAG.EQ.2) THEN
         WRITE (REC,FMT=99997) INAG, NAG
         CALL X04BAF(UNIT1,REC)
      ELSE IF (XNAG.EQ.3) THEN
         WRITE (REC,FMT=99996) INAG, NAG
         CALL X04BAF(UNIT1,REC)
      ELSE IF (XNAG.EQ.4) THEN
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         REC = ' NAG ROUTINE '//NAG
         CALL X04BAF(UNIT1,REC)
         REC = ' ------------------'
         CALL X04BAF(UNIT1,REC)
         CALL X04BAF(UNIT1,' ')
      ELSE
         CALL X04BAF(UNIT1,' ')
         REC = '           * NAG ROUTINE '//NAG//' *'
         CALL X04BAF(UNIT1,REC)
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90HAF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (' ',I3,'.  ',A6,'  -  TEST SUCCESSFULLY PASSED')
99998 FORMAT (' ',I3,'.  ',A6,'  -  **  TEST FAILED (SEE DETAILED OUTP',
     *       'UT BELOW)  **')
99997 FORMAT (' ',I3,'.  ',A6,'  -  **  CANNOT BE TESTED :  CHECK ITS ',
     *       'REQUIREMENTS  **')
99996 FORMAT (' ',I3,'.  ',A6,'  -  **  NOT TESTED :  ROUTINE HAS ALRE',
     *       'ADY PASSED THE TEST  **')
      END
