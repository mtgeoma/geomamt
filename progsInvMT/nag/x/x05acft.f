      INTEGER FUNCTION X05ACF(CTIME1,CTIME2)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Compares the date/time strings CTIME1 and CTIME2.
C     Returns -1 if CTIME1 is earlier than CTIME2,
C              0 if CTIME1 is the same as CTIME2,
C          and 1 if CTIME1 is later than CTIME2.
C
C     .. Scalar Arguments ..
      CHARACTER*(*)           CTIME1, CTIME2
C     .. Local Scalars ..
      INTEGER                 I
C     .. Local Arrays ..
      INTEGER                 ITIME1(7), ITIME2(7)
C     .. External Subroutines ..
      EXTERNAL                X05ACZ
C     .. Executable Statements ..
C
C     Convert each time string into integer array format.
C
      CALL X05ACZ(CTIME1,ITIME1)
      CALL X05ACZ(CTIME2,ITIME2)
C
C     Compare the integer array format times.
C
      I = 1
   20 IF (ITIME1(I).EQ.ITIME2(I) .AND. I.LT.7) THEN
         I = I + 1
         GO TO 20
      END IF
      IF (ITIME1(I).LT.ITIME2(I)) THEN
         X05ACF = -1
      ELSE IF (ITIME1(I).EQ.ITIME2(I)) THEN
         X05ACF = 0
      ELSE
         X05ACF = 1
      END IF
      RETURN
      END
