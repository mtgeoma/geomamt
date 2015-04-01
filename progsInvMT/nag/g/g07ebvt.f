      SUBROUTINE G07EBV(AM,N,X,M,Y,IWRK,SQ)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     ***  Partition step  ***
C
C     Use AM to partition S0 into 2 groups: those .lt. AM, those
C     .ge. AM IWRK(2*N+I)= how many pairs (X(I)+X(J)) in row I
C     less than AM
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AM
      INTEGER           M, N, SQ
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(M)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      INTEGER           I, J, N2
C     .. Executable Statements ..
C
      N2 = 2*N
      J = M
C
C     In upper right corner
C
      SQ = 0
C
C     I counts rows
C
      DO 40 I = 1, N
         IWRK(N2+I) = 0
   20    CONTINUE
C
C        Have we hit the left edge ?
C
         IF (J.LT.1) THEN
            GO TO 40
C
C           Shall we move left ?
C
         ELSE IF (Y(J)-X(N-I+1).GE.AM) THEN
            J = J - 1
            GO TO 20
         END IF
C
C           We're done in this row
C
         IWRK(N2+I) = J
C
C           SQ = Total number of pairs less than AM
C
         SQ = SQ + IWRK(N2+I)
C
   40 CONTINUE
      RETURN
      END
