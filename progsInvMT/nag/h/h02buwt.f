      SUBROUTINE H02BUW(N,NAMES,ID,J1,J2,JMARK,JFOUND)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     ************************************************************
C     H02BUW  FINDS WHERE, IF ANYWHERE, A CERTAIN NAME APPEARS
C     WITHIN A GIVEN LIST OF NAMES.
C     THE CERTAIN NAME IS INPUT IN ARRAY  ID IN A8 FORMAT.
C     THE LIST OF NAMES IS CONTAINED IN ARRAY  NAMES  IN THE SAME
C     FORMAT.
C
C     J1 AND J2 SPECIFY THE LIMITS OF THE SEARCH.
C     ON INPUT,  JMARK SPECIFIES WHERE THE SEARCH SHOULD START.
C     ON OUTPUT, JMARK WILL INDICATE WHERE THE NAME WAS FOUND,
C     OR WILL POINT TO POSITION J1 IF THE NAME WAS NOT FOUND.
C     ON OUTPUT, JFOUND WILL ALSO INDICATE WHERE THE NAME WAS FOUND,
C     BUT WILL BE SET TO ZERO IF THE NAME WAS NOT FOUND.
C
C     ****************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           J1, J2, JFOUND, JMARK, N
      CHARACTER*8       ID
C     .. Array Arguments ..
      CHARACTER*8       NAMES(N)
C     .. Local Scalars ..
      INTEGER           J
C     .. Executable Statements ..
      DO 20 J = JMARK, J2
         IF (ID.EQ.NAMES(J)) GO TO 60
   20 CONTINUE
C
      DO 40 J = J1, JMARK
         IF (ID.EQ.NAMES(J)) GO TO 60
   40 CONTINUE
C
C     *******************************************************
C     Not found
C
      JFOUND = 0
      JMARK = J1
      RETURN
C
C     *******************************************************
C     Found
C
   60 JFOUND = J
      JMARK = J
C
      RETURN
      END
