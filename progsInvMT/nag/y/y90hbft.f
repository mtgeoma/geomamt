      SUBROUTINE Y90HBF(PROG,TPROG,XPROG)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==========================================
C         *  Y90HBF :  Print test program headers  *
C         ==========================================
C
C
C     Print Test program related messages
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           XPROG
      LOGICAL           TPROG
      CHARACTER*6       PROG
C     .. Local Scalars ..
      INTEGER           UNIT1
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Print common header
C
C-----------------------------------------------------------------------
      CALL X04AAF(0,UNIT1)
C
      IF (XPROG.LE.0) THEN
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         REC = ' '//PROG//' STRINGENT TEST PROGRAM RESULTS'
         CALL X04BAF(UNIT1,REC)
         CALL X04BAF(UNIT1,' ')
         REC = ' '
         REC(22:) = '========================='
         CALL X04BAF(UNIT1,REC)
         REC = ' '
         REC(22:) = '*  TEST PROGRAM '//PROG//'  *'
         CALL X04BAF(UNIT1,REC)
         REC = ' '
         REC(22:) = '========================='
         CALL X04BAF(UNIT1,REC)
      ELSE IF (XPROG.EQ.1) THEN
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         REC = ' '
         REC(12:) = 'TEST PROGRAM RESULTS SUMMARY'
         CALL X04BAF(UNIT1,REC)
         REC = ' '
         REC(12:) = '----------------------------'
         CALL X04BAF(UNIT1,REC)
      ELSE IF (XPROG.EQ.2) THEN
         REC(1:) = ' **  TEST PROGRAM  HAS ALREADY BEEN TESTED '//
     *             'SUCCESSFULLY  **'
         CALL X04BAF(UNIT1,REC)
      ELSE IF (XPROG.EQ.3) THEN
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         REC = ' '
         REC(12:) = 'DETAILED DIAGNOSTIC OUTPUT'
         CALL X04BAF(UNIT1,REC)
         REC = ' '
         REC(12:) = '--------------------------'
         CALL X04BAF(UNIT1,REC)
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
      ELSE IF (XPROG.EQ.4) THEN
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         REC = ' NAG ROUTINES INDIVIDUALLY TESTED'
         CALL X04BAF(UNIT1,REC)
         CALL X04BAF(UNIT1,' ')
      ELSE
         CALL X04BAF(UNIT1,' ')
         CALL X04BAF(UNIT1,' ')
         REC = ' JOINT TASKS TESTED'
         CALL X04BAF(UNIT1,REC)
         CALL X04BAF(UNIT1,' ')
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90HBF
C
C-----------------------------------------------------------------------
      RETURN
      END
