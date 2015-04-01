      SUBROUTINE X04CBZ(STRING,START,FINISH)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Returns START as first non blank position of STRING, and FINISH
C     as the last non-blank.
C     .. Scalar Arguments ..
      INTEGER           FINISH, START
      CHARACTER*(*)     STRING
C     .. Local Scalars ..
      INTEGER           L
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (STRING.EQ.' ') THEN
         START = 0
         FINISH = 0
      ELSE
         L = LEN(STRING)
         START = 1
   20    IF (STRING(START:START).EQ.' ' .AND. START.LT.L) THEN
            START = START + 1
            GO TO 20
         END IF
         FINISH = L
   40    IF (STRING(FINISH:FINISH).EQ.' ' .AND. FINISH.GT.1) THEN
            FINISH = FINISH - 1
            GO TO 40
         END IF
      END IF
      RETURN
      END
