      SUBROUTINE G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,POSN,
     *                  CHAR)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Place the character  CHAR at position  POSN  in output line PLOT.
C     If  POSN = 0, place  CHAR in the next available position in  PLOT.
C     MAXPTR is to be initialized to  1, and  it is reset.
C
C     .. Scalar Arguments ..
      INTEGER           CHAR, LDP, MAXPTR, OUTPTR, PMAX, POSN, PSUB
C     .. Array Arguments ..
      CHARACTER         CHARS(*), PLOT(LDP,*)
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      IF (POSN.NE.0) OUTPTR = POSN
      OUTPTR = MIN(OUTPTR,PMAX)
      PLOT(PSUB,OUTPTR) = CHARS(CHAR)
      MAXPTR = MAX(MAXPTR,OUTPTR)
      OUTPTR = OUTPTR + 1
      RETURN
      END
