      SUBROUTINE F11BAZ(ACTION,IDATA,RDATA,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11BAZ - Auxiliary routine for F11GAF, F11BAF, F11BF,
C              F11BBF, F11GCF, F11BCF
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           ACTION, INFO
C     .. Array Arguments ..
      DOUBLE PRECISION  RDATA(20)
      INTEGER           IDATA(20)
C     .. Local Arrays ..
      DOUBLE PRECISION  RDATAX(20)
      INTEGER           IDATAX(20)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06DFF
C     .. Save statement ..
      SAVE              IDATAX, RDATAX
C     .. Data statements ..
      DATA              IDATAX, RDATAX/20*0, 20*ZERO/
C     .. Executable Statements ..
      INFO = 0
C
C     Call from F11GAF, F11BAF
C
      IF (ACTION.EQ.1) THEN
         IDATA(1) = 1
         CALL DCOPY(20,RDATA,1,RDATAX,1)
         CALL F06DFF(20,IDATA,1,IDATAX,1)
C
C     First call from F11GBF, F11BBF
C
      ELSE IF (ACTION.EQ.2) THEN
         IF (IDATAX(1).EQ.1) THEN
            IDATAX(1) = 2
            CALL DCOPY(20,RDATAX,1,RDATA,1)
            CALL F06DFF(20,IDATAX,1,IDATA,1)
         ELSE
            INFO = 1
         END IF
C
C     Subsequent call from F11GBF, F11BBF
C
      ELSE IF (ACTION.EQ.3) THEN
         IF (IDATA(12).EQ.4) IDATA(1) = 3
         CALL DCOPY(20,RDATA,1,RDATAX,1)
         CALL F06DFF(20,IDATA,1,IDATAX,1)
C
C     Call from F11GCF, F11BCF
C
      ELSE IF (ACTION.EQ.4) THEN
         IF ((IDATAX(1).EQ.3) .OR. (IDATAX(12).EQ.3)) THEN
            CALL DCOPY(20,RDATAX,1,RDATA,1)
            CALL F06DFF(20,IDATAX,1,IDATA,1)
         ELSE
            INFO = 1
         END IF
      END IF
C
C     End of subroutine F11BAZ
C
      RETURN
      END
