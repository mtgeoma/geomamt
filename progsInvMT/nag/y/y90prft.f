      SUBROUTINE Y90PRF(LINE,NVECR,VECR,IVECR)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ===================================
C         *  Y90PRF :  Print a Real Vector  *
C         ===================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECR, NVECR
      CHARACTER*(*)     LINE
C     .. Array Arguments ..
      DOUBLE PRECISION  VECR(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DUMMY
      INTEGER           I, IREC, N1, N2, NREC, UNIT1
      CHARACTER*80      REC
C     .. External Functions ..
      INTEGER           X02BEF
      EXTERNAL          X02BEF
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
      REC = ' '
      REC(5:) = LINE(1:)
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
C
      IF (X02BEF(DUMMY).LE.10) THEN
         NREC = (NVECR+2)/3
         N2 = 0
         DO 20 IREC = 1, NREC
            N1 = N2 + 1
            N2 = MIN(N2+3,NVECR)
            WRITE (REC,FMT=99999) (I,VECR((I-1)*IVECR+1),I=N1,N2)
            IF (N2-N1.NE.2) REC(25+(N2-N1)*28:) = ' '
            CALL X04BAF(UNIT1,REC)
   20    CONTINUE
C
      ELSE
         NREC = (NVECR+1)/2
         N2 = 0
         DO 40 IREC = 1, NREC
            N1 = N2 + 1
            N2 = MIN(N2+2,NVECR)
            WRITE (REC,FMT=99998) (I,VECR((I-1)*IVECR+1),I=N1,N2)
            IF (N1.EQ.N2) REC(36:) = ' '
            CALL X04BAF(UNIT1,REC)
   40    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90PRF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,1P,2('(',I3,')  ',D13.6,8X),'(',I3,')  ',D13.6)
99998 FORMAT (4X,1P,2('(',I3,')  ',D19.12,8X))
      END
