      SUBROUTINE Y07ACF(NOUT)
C     .. Scalar Arguments ..
      INTEGER           NOUT
C     .. Local Scalars ..
      CHARACTER*5       REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Executable Statements ..
      WRITE (REC,FMT=99999)
      CALL X04BAF(NOUT,REC)
      RETURN
C
99999 FORMAT (' .+-.')
      END
