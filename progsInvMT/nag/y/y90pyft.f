      SUBROUTINE Y90PYF(LINE,VNAME1,VNAME2,NVECX,VECX1,IVECX1,VECX2,
     *                  IVECX2)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =========================================
C         *  Y90PYF :  Print Two Integer Vectors  *
C         =========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECX1, IVECX2, NVECX
      CHARACTER*(*)     LINE, VNAME1, VNAME2
C     .. Array Arguments ..
      CHARACTER*(*)     VECX1(*), VECX2(*)
C     .. Local Scalars ..
      INTEGER           I, J1, J2, L1, L2, UNIT1
      CHARACTER*10      NAME1, NAME2
      CHARACTER*30      ITEM1, ITEM2
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
      REC(5:) = LINE
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
      NAME1(1:) = VNAME1
      NAME2(1:) = VNAME2
      REC = ' '
      REC(17:) = NAME1
      REC(53:) = NAME2
      CALL X04BAF(UNIT1,REC)
C
      L1 = LEN(VECX1(1))
      L2 = LEN(VECX2(1))
      DO 20 I = 1, NVECX
         J1 = (I-1)*IVECX1 + 1
         J2 = (I-1)*IVECX2 + 1
         ITEM1 = '"'//VECX1(J1)//'"'
         IF (L1.GT.28) ITEM1(27:) = '..."'
         ITEM2 = '"'//VECX2(J2)//'"'
         IF (L2.GT.28) ITEM2(27:) = '..."'
         WRITE (REC,FMT=99999) I, ITEM1, ITEM2
         CALL X04BAF(UNIT1,REC)
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PYF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,'(',I3,')',3X,A30,6X,A30)
      END
