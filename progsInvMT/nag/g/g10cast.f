      SUBROUTINE G10CAS(Y,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     compute repeated running medians of 3
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      LOGICAL           CHANGE
C     .. External Subroutines ..
      EXTERNAL          G10CAR, G10CAW
C     .. Executable Statements ..
C
   20 CHANGE = .FALSE.
      CALL G10CAW(Y,N,CHANGE)
      IF (CHANGE) GO TO 20
      CALL G10CAR(Y,N)
      RETURN
      END
