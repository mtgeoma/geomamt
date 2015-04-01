      SUBROUTINE X05AAF(ITIM)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Returns the current date and time in the integer array ITIM.
C     On exit, ITIM should contain the following values:
C       ITIM(1) : the current year.
C       ITIM(2) : the current month, in the range 1 - 12.
C       ITIM(3) : the current day, in the range 1 - 31.
C       ITIM(4) : the current hour, in the range 0 - 23.
C       ITIM(5) : the current minute, in the range 0 - 59.
C       ITIM(6) : the current second, in the range 0 - 59.
C       ITIM(7) : the current millisecond, in the range 0 - 999.
C
C     This routine is machine dependent. It must be modified by the
C     implementor to return the values described above. If it is not
C     possible on a particular machine to return a correct value for
C     any or all of the above elements, those elements should be set
C     to zero.
C
C     .. Array Arguments ..
      INTEGER           ITIM(7)
C     .. Local Arrays ..
      integer           iarray(3)
C     .. External Subroutines ..
C     EXTERNAL          idate, itime
C     .. Executable Statements ..
C
      call idate(iarray)
      ITIM(1) = iarray(3)
      ITIM(2) = iarray(2)
      ITIM(3) = iarray(1)
      call itime(iarray)
      ITIM(4) = iarray(1)
      ITIM(5) = iarray(2)
      ITIM(6) = iarray(3)
      ITIM(7) = 0
C
      RETURN
      END
