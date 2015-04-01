      SUBROUTINE Y07AEF(NIN,NOUT)
C     .. Scalar Arguments ..
      INTEGER           NIN, NOUT
C     .. Local Scalars ..
      CHARACTER*80      TITLE
      CHARACTER*81      REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Executable Statements ..
      READ (NIN,FMT=99999) TITLE
      WRITE (REC,FMT=99998) TITLE
      CALL X04BAF(NOUT,REC)
      RETURN
C
99999 FORMAT (A80)
99998 FORMAT (1X,A80)
      END
