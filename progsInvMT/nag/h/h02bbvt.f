      SUBROUTINE H02BBV(NXANOD,JM1,J,ZLBAR,AOPTVL,NODTST,NONOD)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     SUBROUTINE CHOOSES AN ACTIVE NODE TO ANALYSE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ZLBAR
      INTEGER           J, JM1, NODTST, NONOD
C     .. Array Arguments ..
      DOUBLE PRECISION  AOPTVL(*)
      INTEGER           NXANOD(0:*)
C     .. Executable Statements ..
      JM1 = 0
      IF (NODTST.GT.NONOD/2) THEN
         J = NXANOD(0)
      ELSE
   20    CONTINUE
         J = NXANOD(JM1)
         IF (AOPTVL(J).EQ.ZLBAR) RETURN
         JM1 = J
         GO TO 20
      END IF
      RETURN
      END
