      SUBROUTINE Y90PLF(LINE,NVECC,VECC,IVECC)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ======================================
C         *  Y90PLF :  Print a Complex Vector  *
C         ======================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECC, NVECC
      CHARACTER*(*)     LINE
C     .. Array Arguments ..
      COMPLEX*16        VECC(*)
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
      INTRINSIC         DIMAG, MIN, DBLE
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
         NREC = (NVECC+1)/2
         N2 = 0
         DO 20 IREC = 1, NREC
            N1 = N2 + 1
            N2 = MIN(N2+2,NVECC)
            WRITE (REC,FMT=99999) (I,VECC((I-1)*IVECC+1),I=N1,N2)
            IF (N2-N1.LE.0) REC(40:) = ' '
            CALL X04BAF(UNIT1,REC)
   20    CONTINUE
C
      ELSE
         NREC = (NVECC+1)/2
         N2 = 0
         DO 40 IREC = 1, NREC
            N1 = N2 + 1
            N2 = MIN(N2+2,NVECC)
            WRITE (REC,FMT=99998) (I,DBLE(VECC((I-1)*IVECC+1)),I=N1,N2)
            IF (N1.EQ.N2) REC(36:65) = ' '
            REC(67:) = 'REAL'
            CALL X04BAF(UNIT1,REC)
            WRITE (REC,FMT=99997) (DIMAG(VECC((I-1)*IVECC+1)),I=N1,N2)
            IF (N1.EQ.N2) REC(36:65) = ' '
            REC(67:) = 'IMAG'
            CALL X04BAF(UNIT1,REC)
   40    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90PLF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,1P,'(',I3,')  ',D13.6,2X,D13.6,6X,'(',I3,')  ',D13.6,
     *       2X,D13.6)
99998 FORMAT (4X,1P,'(',I3,')  ',D19.12,8X,'(',I3,')  ',D19.12)
99997 FORMAT (4X,1P,7X,D19.12,15X,D19.12)
      END
