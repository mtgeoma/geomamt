      SUBROUTINE D02QWZ(ERSTAT,WHATDO)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     TEST THAT SETUP ROUTINE HAS BEEN CALLED
C
C     WHATDO = 1     SETS THE VARIABLE ERSTAT
C            = 0     ASKS WHAT THE LAST VALUE OF ERSTAT WAS
C
C     .. Parameters ..
      INTEGER           SET
      PARAMETER         (SET=1)
C     .. Scalar Arguments ..
      INTEGER           ERSTAT, WHATDO
C     .. Local Scalars ..
      INTEGER           ISAVE
C     .. Save statement ..
      SAVE              ISAVE
C     .. Data statements ..
      DATA              ISAVE/0/
C     .. Executable Statements ..
      IF (WHATDO.EQ.SET) THEN
         ISAVE = ERSTAT
      ELSE
         ERSTAT = ISAVE
      END IF
      RETURN
      END
