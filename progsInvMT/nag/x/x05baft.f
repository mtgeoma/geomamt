      DOUBLE PRECISION FUNCTION X05BAF()
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Returns the amount of CPU time used since some previous time,
C     in seconds. The previous time is arbitrary, but may be the time
C     that the job or the program started, for example. The difference
C     between two separate calls of this routine is then (approximately)
C     the CPU time used between the calls.
C
C     This routine is machine dependent. It must be modified by the
C     implementor to return the value described above.
C
      real etime, tarray(2), secnd
      external etime
      secnd = etime(tarray)
      X05BAF = tarray(1)
      RETURN
      END
