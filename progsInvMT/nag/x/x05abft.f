      CHARACTER*30 FUNCTION X05ABF(ITIME)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Converts from the seven-integer format date/time in ITIME, as
C     returned by subroutine X05AAF, into a character string.
C     The character string has the format, for example,
C     'Thu 20th Apr 1989 16:37:57.320' .
C
C     .. Parameters ..
      INTEGER                      NDAYS, NMNTHS
      PARAMETER                    (NDAYS=7,NMNTHS=12)
C     .. Array Arguments ..
      INTEGER                      ITIME(7)
C     .. Local Scalars ..
      INTEGER                      CENTS, DAY, DOMLET, HOUR, LEAPS,
     *                             MILLI, MINUTE, MONTH, SECOND, THEDAY,
     *                             XLEAPS, YEAR
      LOGICAL                      ISLEAP
      CHARACTER*4                  CHDAY
      CHARACTER*17                 THETIM
C     .. Local Arrays ..
      INTEGER                      MONNOS(NMNTHS), MTHDYS(NMNTHS)
      CHARACTER*3                  DAYS(NDAYS), MONTHS(NMNTHS)
C     .. External Subroutines ..
      EXTERNAL                     X05ABZ
C     .. Intrinsic Functions ..
      INTRINSIC                    MOD
C     .. Data statements ..
      DATA                         MONTHS/'Jan', 'Feb', 'Mar', 'Apr',
     *                             'May', 'Jun', 'Jul', 'Aug', 'Sep',
     *                             'Oct', 'Nov', 'Dec'/
      DATA                         DAYS/'Sun', 'Mon', 'Tue', 'Wed',
     *                             'Thu', 'Fri', 'Sat'/
      DATA                         MONNOS/0, 3, 3, 6, 1, 4, 6, 2, 5, 0,
     *                             3, 5/
      DATA                         MTHDYS/31, 29, 31, 30, 31, 30, 31,
     *                             31, 30, 31, 30, 31/
C     .. Executable Statements ..
      YEAR = ITIME(1)
      MONTH = ITIME(2)
      DAY = ITIME(3)
      HOUR = ITIME(4)
      MINUTE = ITIME(5)
      SECOND = ITIME(6)
      MILLI = ITIME(7)
      LEAPS = YEAR/4
      CENTS = YEAR/100
      XLEAPS = CENTS/4
C
C     Set ISLEAP true if the current year is a leap year.
C
      ISLEAP = MOD(YEAR,4) .EQ. 0
      IF (ISLEAP) THEN
         IF (YEAR.EQ.CENTS*100) THEN
            ISLEAP = MOD(CENTS,4) .EQ. 0
         END IF
      END IF
C
C     Check for illegal date.
C
      IF (MONTH.LT.1 .OR. MONTH.GT.12) THEN
         X05ABF = '** Illegal date **'
      ELSE IF (YEAR.LT.1 .OR. (DAY.LT.1 .OR. DAY.GT.MTHDYS(MONTH))
     *         .OR. (MONTH.EQ.2 .AND. DAY.EQ.29 .AND. .NOT. ISLEAP)
     *         .OR. (HOUR.LT.0 .OR. HOUR.GT.23)
     *         .OR. (MINUTE.LT.0 .OR. MINUTE.GT.59)
     *         .OR. (SECOND.LT.0 .OR. SECOND.GT.59)
     *         .OR. (MILLI.LT.0 .OR. MILLI.GT.999)) THEN
         X05ABF = '** Illegal date **'
      ELSE
C
C        Compute DOMLET, the dominical letter for the year YEAR.
C        1=A, 2=B, 3=C, 4=D, 5=E, 6=F, 7=G.
C
         DOMLET = YEAR + LEAPS - CENTS + XLEAPS
         IF (ISLEAP .AND. MONTH.LE.2) DOMLET = DOMLET - 1
         DOMLET = 8 - MOD(DOMLET,7)
         IF (DOMLET.EQ.8) DOMLET = 1
C
C        Array MONNOS contains offsets for each month, to vector
C        into the correct day-of-the-week table.
C
         THEDAY = 8 - DOMLET + MONNOS(MONTH) + DAY
         THEDAY = MOD(THEDAY,7)
         IF (THEDAY.EQ.0) THEDAY = 7
C
C        THEDAY is the day of the week: 1 = Sunday, 2 = Monday, etc.
C        See Encyclopaedia Brittanica under 'Calendar, perpetual' for
C        the table from which this formula was derived.
C
         IF (DAY.EQ.1 .OR. DAY.EQ.21 .OR. DAY.EQ.31) THEN
            CHDAY = '  st'
         ELSE IF (DAY.EQ.2 .OR. DAY.EQ.22) THEN
            CHDAY = '  nd'
         ELSE IF (DAY.EQ.3 .OR. DAY.EQ.23) THEN
            CHDAY = '  rd'
         ELSE
            CHDAY = '  th'
         END IF
C
         IF (DAY.GT.9) THEN
            CALL X05ABZ(CHDAY(1:2),DAY)
         ELSE
            CALL X05ABZ(CHDAY(2:2),DAY)
         END IF
C
         THETIM = '       :  :  .   '
         CALL X05ABZ(THETIM(1:4),YEAR)
         CALL X05ABZ(THETIM(6:7),HOUR)
         CALL X05ABZ(THETIM(9:10),MINUTE)
         CALL X05ABZ(THETIM(12:13),SECOND)
         CALL X05ABZ(THETIM(15:17),MILLI)
C
         X05ABF = DAYS(THEDAY)//' '//CHDAY//' '//MONTHS(MONTH)
     *            //' '//THETIM
      END IF
C
      RETURN
      END
