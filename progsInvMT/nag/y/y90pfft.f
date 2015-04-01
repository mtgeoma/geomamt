      SUBROUTINE Y90PFF(LINE,CVAR)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90PFF :  Print a Complex Variable  *
C         =======================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      COMPLEX*16        CVAR
      CHARACTER*(*)     LINE
C     .. Local Scalars ..
      DOUBLE PRECISION  DUMMY
      INTEGER           UNIT1
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
      REC(50:) = ' :  '
      IF (X02BEF(DUMMY).LE.10) THEN
         WRITE (REC(54:),FMT=99999) CVAR
         CALL X04BAF(UNIT1,REC)
      ELSE
         WRITE (REC(54:),FMT=99998) DBLE(CVAR)
         CALL X04BAF(UNIT1,REC)
         WRITE (REC,FMT=99997) DIMAG(CVAR)
         CALL X04BAF(UNIT1,REC)
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90PFF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (1P,D13.6,1X,D13.6)
99998 FORMAT (1P,D19.12,' REAL')
99997 FORMAT (1P,53X,D19.12,' IMAG')
      END
