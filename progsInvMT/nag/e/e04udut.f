      SUBROUTINE E04UDU(STRING)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     PURPOSE:  This subroutine changes all lower case letters in the
C               character string to upper case.
C
C     METHOD:   Each character in STRING is treated in turn.  The
C               intrinsic function INDEX effectively allows a table
C               lookup, with the local strings LOW and UPP acting as
C               two tables. This method avoids the use of CHAR and
C               ICHAR, which appear be different on ASCII and EBCDIC
C               machines.
C
C     ARGUMENTS
C     ARG       DIM     TYPE I/O/S DESCRIPTION
C     STRING       *       C   I/O   Character string possibly
C                                   containing some lower-case
C                                   letters  on input; strictly
C                                   upper-case letters on output
C                                   with no change to any
C                                   non-alphabetic characters.
C
C     EXTERNAL REFERENCES:
C     LEN    - Returns the declared length of a CHARACTER variable.
C     INDEX  - Returns the position of second string within first.
C
C     ENVIRONMENT:  ANSI FORTRAN 77
C
C     DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C     06/28/83   CLH    Initial design.
C     01/03/84   RAK    Eliminated NCHAR input.
C     06/14/84   RAK    Used integer PARAMETERs in comparison.
C     04/21/85   RAK    Eliminated DO/END DO in favor of standard code.
C     09/10/85   MAS    Eliminated CHAR,ICHAR in favor of LOW, UPP,
C                       INDEX.
C
C     AUTHOR: Charles Hooper, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      CHARACTER*(*)     STRING
C     .. Local Scalars ..
      INTEGER           I, J
      CHARACTER*1       C
      CHARACTER*26      LOW, UPP
C     .. Intrinsic Functions ..
      INTRINSIC         INDEX, LEN, LGE, LLE
C     .. Data statements ..
      DATA              LOW/'abcdefghijklmnopqrstuvwxyz'/,
     *                  UPP/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C     .. Executable Statements ..
C
      DO 20 J = 1, LEN(STRING)
         C = STRING(J:J)
         IF (LGE(C,'a') .AND. LLE(C,'z')) THEN
C           IF (C.GE.'a' .AND. C.LE.'z') THEN
            I = INDEX(LOW,C)
            IF (I.GT.0) STRING(J:J) = UPP(I:I)
         END IF
   20 CONTINUE
      RETURN
C
C     End of  E04UDU. (OPUPPR)
C
      END
