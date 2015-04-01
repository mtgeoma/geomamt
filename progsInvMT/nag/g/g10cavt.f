      SUBROUTINE G10CAV(Y,N,ENDSAV,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     smooth by running medians of 4
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ENDSAV
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(4), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ENDM1
C     .. Local Arrays ..
      INTEGER           IWORK(4)
C     .. External Subroutines ..
      EXTERNAL          G10CAN
C     .. Executable Statements ..
C
C     even length medians offset the output sequence to the upper end,
C     ENDSAV holds Y(N), Y(1) is unchanged
C
      ENDSAV = Y(N)
      ENDM1 = Y(N-1)
      CALL G10CAN(Y,N,4,WORK,IWORK)
C
      Y(2) = (Y(1)+Y(2))/2.0D0
      Y(N) = (ENDM1+ENDSAV)/2.0D0
      RETURN
      END
