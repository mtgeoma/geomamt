      SUBROUTINE Y90PUF(LINE,VNAME1,VNAME2,NVECL,VECL1,IVECL1,VECL2,
     *                  IVECL2)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =========================================
C         *  Y90PUF :  Print Two Logical Vectors  *
C         =========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECL1, IVECL2, NVECL
      CHARACTER*(*)     LINE, VNAME1, VNAME2
C     .. Array Arguments ..
      LOGICAL           VECL1(*), VECL2(*)
C     .. Local Scalars ..
      INTEGER           I, UNIT1
      CHARACTER*8       VALUE1, VALUE2
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
      NAME1(1:) = VNAME1
      NAME2(1:) = VNAME2
C
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
      REC = ' '
      REC(22:) = NAME1
      REC(45:) = NAME2
      CALL X04BAF(UNIT1,REC)
      DO 20 I = 1, NVECL
         IF (VECL1((I-1)*IVECL1+1)) THEN
            VALUE1 = ' .TRUE.'
         ELSE
            VALUE1 = ' .FALSE.'
         END IF
         IF (VECL2((I-1)*IVECL2+1)) THEN
            VALUE2 = ' .TRUE.'
         ELSE
            VALUE2 = ' .FALSE.'
         END IF
         WRITE (REC,FMT=99999) I, VALUE1, VALUE2
         CALL X04BAF(UNIT1,REC)
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PUF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,'(',I3,')',12X,A8,15X,A8)
      END
