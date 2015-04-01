      SUBROUTINE G10CAY(Y,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Smooth  by 3RSSH, twice
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      LOGICAL           CHANGE
C     .. External Subroutines ..
      EXTERNAL          G10CAQ, G10CAS, G10CAT
C     .. Executable Statements ..
C
      CALL G10CAS(Y,N)
      CHANGE = .FALSE.
      CALL G10CAQ(Y,N,CHANGE)
      IF ( .NOT. CHANGE) GO TO 20
      CALL G10CAS(Y,N)
      CHANGE = .FALSE.
      CALL G10CAQ(Y,N,CHANGE)
      IF (CHANGE) CALL G10CAS(Y,N)
   20 CALL G10CAT(Y,N)
      RETURN
      END
