      SUBROUTINE G10CAW(Y,N,CHANGE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     compute running median of 3 on Y().
C     sets CHANGE  .TRUE.  if any CHANGE is made
C
C     .. Scalar Arguments ..
      INTEGER           N
      LOGICAL           CHANGE
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y1, Y2, Y3
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          G10CAP
C     .. Executable Statements ..
C
      Y2 = Y(1)
      Y3 = Y(2)
      DO 20 I = 2, N - 1
         Y1 = Y2
         Y2 = Y3
         Y3 = Y(I+1)
         CALL G10CAP(Y1,Y2,Y3,Y(I),CHANGE)
   20 CONTINUE
      RETURN
      END
