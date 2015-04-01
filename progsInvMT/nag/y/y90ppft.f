      SUBROUTINE Y90PPF(LINE,VNAME1,VNAME2,NVECI,VECI1,IVECI1,VECI2,
     *                  IVECI2)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =========================================
C         *  Y90PPF :  Print Two Integer Vectors  *
C         =========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECI1, IVECI2, NVECI
      CHARACTER*(*)     LINE, VNAME1, VNAME2
C     .. Array Arguments ..
      INTEGER           VECI1(*), VECI2(*)
C     .. Local Scalars ..
      INTEGER           I, J1, J2, UNIT1
      CHARACTER*10      NAME1, NAME2
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
      REC(5:) = LINE
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
      NAME1(1:) = VNAME1
      NAME2(1:) = VNAME2
      REC = ' '
      REC(22:) = NAME1
      REC(45:) = NAME2
      CALL X04BAF(UNIT1,REC)
C
      DO 20 I = 1, NVECI
         J1 = (I-1)*IVECI1 + 1
         J2 = (I-1)*IVECI2 + 1
         WRITE (REC,FMT=99999) I, VECI1(J1), VECI2(J2)
         CALL X04BAF(UNIT1,REC)
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PPF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,'(',I3,')',12X,I8,15X,I8)
      END
