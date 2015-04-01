      SUBROUTINE X04BAY(NOUT,NREC,REC)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAY outputs NREC records on device NOUT, by calling X04BAF.
C     If NREC is 0 then no records are output.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT, NREC
C     .. Array Arguments ..
      CHARACTER*(*)     REC(*)
C     .. Local Scalars ..
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Executable Statements ..
      DO 20 I = 1, NREC
         CALL X04BAF(NOUT,REC(I))
   20 CONTINUE
      RETURN
C
C     End of X04BAY.
C
      END
