      SUBROUTINE D02NXF(ICALL,LIWREQ,LIWUSD,LRWREQ,LRWUSD,NLU,NNZ,NGP,
     *                  ISPLIT,IGROW,LBLOCK,NBLOCK,INFORM)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-712 (DEC 1989).
C     .. Scalar Arguments ..
      INTEGER           ICALL, IGROW, ISPLIT, LIWREQ, LIWUSD, LRWREQ,
     *                  LRWUSD, NBLOCK, NGP, NLU, NNZ
      LOGICAL           LBLOCK
C     .. Array Arguments ..
      INTEGER           INFORM(23)
C     .. Executable Statements ..
      LIWREQ = INFORM(10)
      LIWUSD = INFORM(11)
      LRWREQ = INFORM(12)
      LRWUSD = INFORM(13)
      IF (ICALL.EQ.1) RETURN
      NLU = INFORM(14)
      NNZ = INFORM(15)
      NGP = INFORM(16)
      ISPLIT = INFORM(17)
      IGROW = INFORM(18)
      NBLOCK = 1
      IF (LBLOCK) NBLOCK = INFORM(19)
      RETURN
      END
