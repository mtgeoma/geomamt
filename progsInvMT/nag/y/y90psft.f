      SUBROUTINE Y90PSF(LINE,VNAME1,VNAME2,NVECR,VECR1,IVECR1,VECR2,
     *                  IVECR2)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ======================================
C         *  Y90PSF :  Print Two Real Vectors  *
C         ======================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECR1, IVECR2, NVECR
      CHARACTER*(*)     LINE, VNAME1, VNAME2
C     .. Array Arguments ..
      DOUBLE PRECISION  VECR1(*), VECR2(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DUMMY
      INTEGER           I, J1, J2, UNIT1
      CHARACTER*10      NAME1, NAME2
      CHARACTER*80      REC
C     .. External Functions ..
      INTEGER           X02BEF
      EXTERNAL          X02BEF
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
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
C
      REC = ' '
      IF (X02BEF(DUMMY).LE.10) THEN
         REC(22:) = NAME1
         REC(45:) = NAME2
         CALL X04BAF(UNIT1,REC)
         DO 20 I = 1, NVECR
            J1 = (I-1)*IVECR1 + 1
            J2 = (I-1)*IVECR2 + 1
            WRITE (REC,FMT=99999) I, VECR1(J1), VECR2(J2)
            CALL X04BAF(UNIT1,REC)
   20    CONTINUE
      ELSE
         REC(24:) = NAME1
         REC(53:) = NAME2
         CALL X04BAF(UNIT1,REC)
         DO 40 I = 1, NVECR
            J1 = (I-1)*IVECR1 + 1
            J2 = (I-1)*IVECR2 + 1
            WRITE (REC,FMT=99998) I, VECR1(J1), VECR2(J2)
            CALL X04BAF(UNIT1,REC)
   40    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90PSF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,'(',I3,')',1P,2(10X,D13.6))
99998 FORMAT (4X,'(',I3,')',1P,2(10X,D19.12))
      END
