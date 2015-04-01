      LOGICAL FUNCTION E04UDX(STRING)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C***********************************************************************
C     Description and usage:
C
C        A simple(-minded) test for numeric data is implemented by
C        searching an input string for legitimate characters:
C                digits 0 to 9, D, E, -, + and .
C        Insurance is provided by requiring that a numeric string
C        have at least one digit, at most one D, E or .
C        and at most two -s or +s.  Note that a few ambiguities remain:
C
C           (a)  A string might have the form of numeric data but be
C                intended as text.  No general test can hope to detect
C                such cases.
C
C           (b)  There is no check for correctness of the data format.
C                For example a meaningless string such as 'E1.+2-'
C                will be accepted as numeric.
C
C        Despite these weaknesses, the method should work in the
C        majority of cases.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        E04UDX              L      O    Set .TRUE. if STRING appears
C                                        to be numerical data.
C        STRING              C    I      Input data to be tested.
C
C
C     Environment:  ANSI FORTRAN 77.
C
C
C     Notes:
C
C        (1)  It is assumed that STRING is a token extracted by
C             E04UDV, which will have converted any lower-case
C             characters to upper-case.
C
C        (2)  E04UDV pads STRING with blanks, so that a genuine
C             number is of the form  '1234        '.
C             Hence, the scan of STRING stops at the first blank.
C
C        (3)  COMPLEX data with parentheses will not look numeric.
C
C
C     Systems Optimization Laboratory, Stanford University.
C     12 Nov  1985    Initial design and coding, starting from the
C                     routine ALPHA from Informatics General, Inc.
C***********************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER*(*)           STRING
C     .. Local Scalars ..
      INTEGER                 J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS,
     *                        NPOINT
      LOGICAL                 NUMBER
      CHARACTER*1             ATOM
C     .. Intrinsic Functions ..
      INTRINSIC               LEN, LGE, LLE
C     .. Executable Statements ..
      NDIGIT = 0
      NEXP = 0
      NMINUS = 0
      NPLUS = 0
      NPOINT = 0
      NUMBER = .TRUE.
      LENGTH = LEN(STRING)
      J = 0
C
   20 J = J + 1
      ATOM = STRING(J:J)
      IF (LGE(ATOM,'0') .AND. LLE(ATOM,'9')) THEN
C        IF (ATOM.GE.'0' .AND. ATOM.LE.'9') THEN
         NDIGIT = NDIGIT + 1
      ELSE IF (ATOM.EQ.'D' .OR. ATOM.EQ.'E') THEN
         NEXP = NEXP + 1
      ELSE IF (ATOM.EQ.'-') THEN
         NMINUS = NMINUS + 1
      ELSE IF (ATOM.EQ.'+') THEN
         NPLUS = NPLUS + 1
      ELSE IF (ATOM.EQ.'.') THEN
         NPOINT = NPOINT + 1
      ELSE IF (ATOM.EQ.' ') THEN
         J = LENGTH
      ELSE
         NUMBER = .FALSE.
      END IF
C
      IF (NUMBER .AND. J.LT.LENGTH) GO TO 20
C
      E04UDX = NUMBER .AND. NDIGIT .GE. 1 .AND. NEXP .LE. 1 .AND.
     *         NMINUS .LE. 2 .AND. NPLUS .LE. 2 .AND. NPOINT .LE. 1
C
      RETURN
C
C     End of  E04UDX. (OPNUMB)
      END
