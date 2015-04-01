      SUBROUTINE E04UDV(STRING,NUMIN,NUMOUT,LIST)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Description and usage:
C
C       An aid to parsing input data.  The individual 'tokens' in a
C     character string are isolated, converted to uppercase, and stored
C     in an array.  Here, a token is a group of significant, contiguous
C     characters.  The following are NON-significant, and hence may
C     serve as separators:  blanks, horizontal tabs, commas, colons,
C     and equal signs.  See E04UDW for details.  Processing continues
C     until the requested number of tokens have been found or the end
C     of the input string is reached.
C
C
C     Parameters:
C
C     Name    Dimension  Type  I/O/S  Description
C     STRING              C    I      Input string to be analyzed.
C     NUMIN               I    I/O    Number of tokens requested (input)
C     NUMOUT                          and found (output).
C     (NUMIN and NUMOUT were both called NUMBER in the original)
C
C     LIST    NUMIN       C      O    Array of tokens, changed to upper
C                                    case.
C
C
C     External references:
C
C     Name    Description
C     E04UDW  Finds positions of first and last significant characters.
C     E04UDU  Converts a string to uppercase.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C     Notes:
C
C     (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C     16 Jan. 1984    RAK    Initial design and coding.
C     16 Mar. 1984    RAK    Revised header to reflect full list of
C                            separators, repaired faulty WHILE clause
C                            in '10' loop.
C     18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                            at a time, leaving STRING unchanged.
C
C-----------------------------------------------------------------------
C
C     .. Parameters ..
      CHARACTER         BLANK
      PARAMETER         (BLANK=' ')
C     .. Scalar Arguments ..
      INTEGER           NUMIN, NUMOUT
      CHARACTER*(*)     STRING
C     .. Array Arguments ..
      CHARACTER*(*)     LIST(NUMIN)
C     .. Local Scalars ..
      INTEGER           COUNT, FIRST, I, LAST, MARK
C     .. External Subroutines ..
      EXTERNAL          E04UDU, E04UDW
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      FIRST = 1
      LAST = LEN(STRING)
C
      COUNT = 0
   20 CONTINUE
C
C        Get delimiting indices of next token, if any.
C
      CALL E04UDW(STRING,FIRST,LAST,MARK)
      IF (LAST.GT.0) THEN
         COUNT = COUNT + 1
C
C           Pass token to output string array, then change case.
C
         LIST(COUNT) = STRING(FIRST:MARK)
         CALL E04UDU(LIST(COUNT))
         FIRST = MARK + 2
         IF (COUNT.LT.NUMIN) GO TO 20
C
      END IF
C
C
C     Fill the rest of LIST with blanks and set NUMBER for output.
C
      DO 40 I = COUNT + 1, NUMIN
         LIST(I) = BLANK
   40 CONTINUE
C
      NUMOUT = COUNT
C
C
C     Termination.
C     ------------
C
      RETURN
C
C     End of  E04UDV. (OPTOKN)
      END
