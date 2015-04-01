      SUBROUTINE D02LAX(STATE,WHATDO)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     TEST THAT INTEGRATOR ROUTINE HAS BEEN CALLED
C
C     WHATDO = 1     SETS THE VARIABLE STATE
C            = 0     ASKS WHAT THE LAST VALUE OF STATE WAS
C
C     .. Parameters ..
      INTEGER           SET
      PARAMETER         (SET=1)
C     .. Scalar Arguments ..
      INTEGER           STATE, WHATDO
C     .. Local Scalars ..
      INTEGER           ISAVE
C     .. Save statement ..
      SAVE              ISAVE
C     .. Data statements ..
      DATA              ISAVE/0/
C     .. Executable Statements ..
      IF (WHATDO.EQ.SET) THEN
         ISAVE = STATE
      ELSE
         STATE = ISAVE
      END IF
      RETURN
      END
