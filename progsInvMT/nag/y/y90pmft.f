      SUBROUTINE Y90PMF(LINE,VNAME1,VNAME2,NVECC,VECC1,IVECC1,VECC2,
     *                  IVECC2)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =========================================
C         *  Y90PMF :  Print Two Complex Vectors  *
C         =========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECC1, IVECC2, NVECC
      CHARACTER*(*)     LINE, VNAME1, VNAME2
C     .. Array Arguments ..
      COMPLEX*16        VECC1(*), VECC2(*)
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
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, DBLE
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
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
      NAME1(1:) = VNAME1
      NAME2(1:) = VNAME2
C
      REC = ' '
      IF (X02BEF(DUMMY).LE.10) THEN
         REC(26:) = NAME1
         REC(61:) = NAME2
         CALL X04BAF(UNIT1,REC)
         DO 20 I = 1, NVECC
            J1 = (I-1)*IVECC1 + 1
            J2 = (I-1)*IVECC2 + 1
            WRITE (REC,FMT=99999) I, VECC1(J1), VECC2(J2)
            CALL X04BAF(UNIT1,REC)
   20    CONTINUE
C
      ELSE
         REC(24:) = NAME1
         REC(53:) = NAME2
         CALL X04BAF(UNIT1,REC)
         DO 40 I = 1, NVECC
            J1 = (I-1)*IVECC1 + 1
            J2 = (I-1)*IVECC2 + 1
            WRITE (REC,FMT=99998) I, DBLE(VECC1(J1)), DBLE(VECC2(J2))
            CALL X04BAF(UNIT1,REC)
            WRITE (REC,FMT=99997) DIMAG(VECC1(J1)), DIMAG(VECC2(J2))
            CALL X04BAF(UNIT1,REC)
   40    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90PMF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,'(',I3,')',1P,2(7X,D13.6,2X,D13.6))
99998 FORMAT (4X,'(',I3,')',1P,2(10X,D19.12),'  REAL')
99997 FORMAT (9X,1P,2(10X,D19.12),'  IMAG')
      END
