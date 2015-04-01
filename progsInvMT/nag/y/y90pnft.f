      SUBROUTINE Y90PNF(LINE,NVECI,VECI,IVECI)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90PNF :  Print an Integer Vector  *
C         =======================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECI, NVECI
      CHARACTER*(*)     LINE
C     .. Array Arguments ..
      INTEGER           VECI(*)
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
      REC = ' '
      REC(5:) = LINE(1:)
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
      NREC = (NVECI+3)/4
      N2 = 0
      DO 20 IREC = 1, NREC
         N1 = N2 + 1
         N2 = MIN(N2+4,NVECI)
         WRITE (REC,FMT=99999) (I,VECI((I-1)*IVECI+1),I=N1,N2)
         IF (N2-N1.NE.3) REC(24+(N2-N1)*19:) = ' '
         CALL X04BAF(UNIT1,REC)
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PNF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,4('(',I3,')',I8,6X))
      END
