      SUBROUTINE Y90HCF(ITASK,TASK,TTASK,WTASK,XTASK,NAG,NNAG)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ====================================================
C         *  Y90HCF :  Print message for Joint Task testing  *
C         ====================================================
C
C
C     Print Joint Task Headers
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           ITASK, NNAG, XTASK
      LOGICAL           TTASK, WTASK
      CHARACTER*6       TASK
C     .. Array Arguments ..
      CHARACTER*6       NAG(NNAG)
C     .. Local Scalars ..
      INTEGER           I, IREC, N1, N2, NREC, UNIT1
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Print
C
C-----------------------------------------------------------------------
      CALL X04AAF(0,UNIT1)
C
      IF (XTASK.LE.1) THEN
         IF (TTASK) THEN
            WRITE (REC,FMT=99999) ITASK, TASK
            CALL X04BAF(UNIT1,REC)
         ELSE
            WRITE (REC,FMT=99998) ITASK, TASK
            CALL X04BAF(UNIT1,REC)
         END IF
         IF ( .NOT. WTASK) THEN
            REC = ' '
            REC(19:) = '++  WARNINGS HAVE BEEN RAISED'
            CALL X04BAF(UNIT1,REC)
         END IF
C
         N1 = MIN(5,NNAG)
         WRITE (REC,FMT=99997) (NAG(I),I=1,N1)
         CALL X04BAF(UNIT1,REC)
         NREC = (NNAG-N1+4)/5
         DO 20 IREC = 1, NREC
            N1 = 5*IREC + 1
            N2 = MIN(NNAG,N1+4)
            WRITE (REC,FMT=99996) (NAG(I),I=N1,N2)
            CALL X04BAF(UNIT1,REC)
   20    CONTINUE
C
      ELSE IF (XTASK.EQ.2) THEN
         WRITE (REC,FMT=99995) ITASK, TASK
         CALL X04BAF(UNIT1,REC)
      ELSE IF (XTASK.EQ.3) THEN
         WRITE (REC,FMT=99994) ITASK, TASK
         CALL X04BAF(UNIT1,REC)
      ELSE IF (XTASK.EQ.4) THEN
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         REC = ' JOINT TASK '//TASK
         CALL X04BAF(UNIT1,REC)
         REC = ' -----------------'
         CALL X04BAF(UNIT1,REC)
      ELSE
         CALL X04BAF(UNIT1,' ')
         REC = '           * JOINT TASK '//TASK//' *'
         CALL X04BAF(UNIT1,REC)
         CALL X04BAF(UNIT1,' ')
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90HCF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (' ',I3,'.  ',A6,'  -  TEST SUCCESSFULLY PASSED')
99998 FORMAT (' ',I3,'.  ',A6,'  -  **  TEST FAILED (SEE DETAILED OUTP',
     *       'UT BELOW)  **')
99997 FORMAT (18X,'NAG ROUTINES :',5(2X,A6))
99996 FORMAT (32X,5(2X,A6))
99995 FORMAT (' ',I3,'.  ',A6,'  -  **  CANNOT BE TESTED :  CHECK ITS ',
     *       'REQUIREMENTS  **')
99994 FORMAT (' ',I3,'.  ',A6,'  -  **  NOT TESTED :  TASK HAS ALREADY',
     *       ' PASSED THE TEST SUCCESSFULLY  **')
      END
