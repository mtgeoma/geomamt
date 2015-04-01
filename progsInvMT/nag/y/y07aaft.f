      SUBROUTINE Y07AAF(NOUT)
C     .. Scalar Arguments ..
      INTEGER           NOUT
C     .. Local Scalars ..
      CHARACTER*5       REC
C     .. External Subroutines ..
      EXTERNAL          A00AAF, X04BAF
C     .. Executable Statements ..
      WRITE (REC,FMT=99999)
      CALL X04BAF(NOUT,REC)
      CALL A00AAF
      RETURN
C
99999 FORMAT (' .++.')
      END
