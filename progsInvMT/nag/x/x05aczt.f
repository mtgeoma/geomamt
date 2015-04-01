      SUBROUTINE X05ACZ(CTIME,ITIME)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-934 (APR 1991).
C
C     Converts from character string format time to integer array
C     format time.
C     Subroutine X05ACZ prefers CTIME to be in the CHARACTER*30
C     format returned by routine X05ABF, but will make an attempt
C     to interpret any string. Any day-name at the start of the
C     string is ignored. Each field should be separated by at
C     least one blank. Thus X05ACZ will regard all the following
C     strings as equivalent:
C
C       'Thu 20th Apr 1989 16:37:57.320'
C       'Thursday 20th Apr 1989 16:37:57.320'
C       'Thursday 20th April 1989 16:37:57.320'
C       '20th April 1989 16:37:57.320'
C
C     .. Parameters ..
      INTEGER           NMNTHS
      PARAMETER         (NMNTHS=12)
C     .. Scalar Arguments ..
      CHARACTER*(*)     CTIME
C     .. Array Arguments ..
      INTEGER           ITIME(7)
C     .. Local Scalars ..
      INTEGER           I, WRDLEN
      CHARACTER*80      STIME, WORD
C     .. Local Arrays ..
      CHARACTER*9       MONTHS(NMNTHS)
C     .. External Functions ..
      INTEGER           X05AAZ
      LOGICAL           X05ACY
      CHARACTER*80      X05ACX
      EXTERNAL          X05AAZ, X05ACY, X05ACX
C     .. Intrinsic Functions ..
      INTRINSIC         LGE, LGT, LLE, LLT
C     .. Data statements ..
      DATA              MONTHS/'JANUARY', 'FEBRUARY', 'MARCH', 'APRIL',
     *                  'MAY', 'JUNE', 'JULY', 'AUGUST', 'SEPTEMBER',
     *                  'OCTOBER', 'NOVEMBER', 'DECEMBER'/
C     .. Executable Statements ..
      STIME = CTIME
      WORD = X05ACX(STIME,WRDLEN)
      IF (WRDLEN.GE.1) THEN
         IF (LLT(WORD(1:1),'0') .OR. LGT(WORD(1:1),'9')) THEN
C
C           Assume that we have an alphabetic day name here; ignore it.
C
            WORD = X05ACX(STIME,WRDLEN)
         END IF
      END IF
C
C     Now WORD contains what we expect to be a day number.
C     Allow the last few characters of WORD to be non-numeric,
C     e.g. WORD could be '12th' or '21st'. On the other hand,
C     don't enforce these non-numeric characters. Thus
C     '12' is equivalent to '12th' and '12buckle-my-shoe'.
C
   20 IF (WRDLEN.GE.1) THEN
         IF (LLT(WORD(WRDLEN:WRDLEN),'0') .OR.
     *       LGT(WORD(WRDLEN:WRDLEN),'9')) THEN
            WRDLEN = WRDLEN - 1
            GO TO 20
         END IF
      END IF
C
      IF (WRDLEN.GE.1) THEN
         ITIME(3) = X05AAZ(WORD(1:WRDLEN))
      ELSE
         ITIME(3) = 0
      END IF
C
      WORD = X05ACX(STIME,WRDLEN)
C
C     Now we expect WORD to contain a month-name, possibly abbreviated.
C     However, we also allow a month number.
C
      IF (WRDLEN.GE.1) THEN
         IF (LGE(WORD(1:1),'0') .AND. LLE(WORD(1:1),'9')) THEN
C
C           The month has been supplied as a number.
C
            ITIME(2) = X05AAZ(WORD(1:WRDLEN))
         ELSE
C
C           Compare with the list of month names until we find a match.
C
            I = 1
   40       IF (X05ACY(WORD(1:WRDLEN),MONTHS(I))) THEN
               ITIME(2) = I
            ELSE IF (I.LT.NMNTHS) THEN
               I = I + 1
               GO TO 40
            ELSE
               ITIME(2) = 0
            END IF
         END IF
      ELSE
         ITIME(2) = 0
      END IF
C
      WORD = X05ACX(STIME,WRDLEN)
C
C     WORD should now contain a year number.
C
      IF (WRDLEN.GE.1) THEN
         ITIME(1) = X05AAZ(WORD(1:WRDLEN))
      ELSE
         ITIME(1) = 0
      END IF
C
      WORD = X05ACX(STIME,WRDLEN)
C
C     Now WORD should contain a time of up to 12 characters,
C     e.g. '11:49:23.327' . We assume missing fields to be zero.
C
      IF (WRDLEN.LT.2) THEN
         ITIME(4) = 0
      ELSE
         ITIME(4) = X05AAZ(WORD(1:2))
      END IF
C
      IF (WRDLEN.LT.5) THEN
         ITIME(5) = 0
      ELSE
         ITIME(5) = X05AAZ(WORD(4:5))
      END IF
C
      IF (WRDLEN.LT.8) THEN
         ITIME(6) = 0
      ELSE
         ITIME(6) = X05AAZ(WORD(7:8))
      END IF
C
      IF (WRDLEN.LT.12) THEN
         ITIME(7) = 0
      ELSE
         ITIME(7) = X05AAZ(WORD(10:12))
      END IF
      RETURN
      END
