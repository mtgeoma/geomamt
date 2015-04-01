      SUBROUTINE G10CAU(Y,N,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     smooth by running medians of 5
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(5), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  YMED1, YMED2
      LOGICAL           CHANGE
C     .. Local Arrays ..
      INTEGER           IWORK(5)
C     .. External Subroutines ..
      EXTERNAL          G10CAN, G10CAP
C     .. Executable Statements ..
C
      CHANGE = .FALSE.
C
      CALL G10CAP(Y(1),Y(2),Y(3),YMED1,CHANGE)
      CALL G10CAP(Y(N),Y(N-1),Y(N-2),YMED2,CHANGE)
      CALL G10CAN(Y,N,5,WORK,IWORK)
      Y(2) = YMED1
      Y(N-1) = YMED2
      RETURN
      END
