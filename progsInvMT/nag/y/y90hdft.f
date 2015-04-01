      SUBROUTINE Y90HDF(ITEST,TEST,WARN,TTEST)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =============================================
C         *  Y90HDF :  Print Individual Test Heading  *
C         =============================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           ITEST
      LOGICAL           TEST, WARN
      CHARACTER*(*)     TTEST
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
      CALL X04BAF(UNIT1,' ')
      IF (TEST) THEN
         IF (WARN) THEN
            WRITE (REC,FMT=99999) ITEST
         ELSE
            WRITE (REC,FMT=99998) ITEST
         END IF
         CALL X04BAF(UNIT1,REC)
      ELSE
         WRITE (REC,FMT=99997) ITEST
         CALL X04BAF(UNIT1,REC)
         IF ( .NOT. WARN) THEN
            REC = ' '
            REC(17:) = '++  WARNING(S) HAVE BEEN RAISED'
            CALL X04BAF(UNIT1,REC)
         END IF
      END IF
      REC = ' '//TTEST
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
C-----------------------------------------------------------------------
C
C     End of Y90HDF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (' TEST NO.',I3,' :  PASSED SUCCESSFULLY')
99998 FORMAT (' TEST NO.',I3,' :  NON-FATAL WARNING(S) RAISED')
99997 FORMAT (' TEST NO.',I3,' :  **  FAILED  **')
      END
