      SUBROUTINE G10CAZ(Y,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Smooth by 4253H
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ENDSAV
      LOGICAL           CHANGE
C     .. Local Arrays ..
      DOUBLE PRECISION  WORK(10)
C     .. External Subroutines ..
      EXTERNAL          G10CAR, G10CAT, G10CAU, G10CAV, G10CAW, G10CAX
C     .. Executable Statements ..
C
      CHANGE = .FALSE.
C
      CALL G10CAV(Y,N,ENDSAV,WORK)
      CALL G10CAX(Y,N,ENDSAV)
      CALL G10CAU(Y,N,WORK)
      CALL G10CAW(Y,N,CHANGE)
      CALL G10CAR(Y,N)
      CALL G10CAT(Y,N)
      RETURN
      END
