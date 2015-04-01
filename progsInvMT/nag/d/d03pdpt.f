      SUBROUTINE D03PDP(ICALLD,IDO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           ICALLD, IDO
C     .. Local Scalars ..
      INTEGER           ISAVE
C     .. Save statement ..
      SAVE              ISAVE
C     .. Data statements ..
      DATA              ISAVE/0/
C     .. Executable Statements ..
      IF (IDO.EQ.1) THEN
         ICALLD = ISAVE
      ELSE
         ISAVE = ICALLD
      END IF
      RETURN
      END
