      SUBROUTINE G07EAV(AM,N,X,IWRK,SQ)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Auxillary routine for G07EAW which uses the exact method for
C     computing the Hodges-Lehmann estimator and the corresponding
C     confidence interval.
C
C     This routine performs the partition step.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AM
      INTEGER           N, SQ
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      INTEGER           I, J, N2
C     .. Executable Statements ..
C
C     Use AM to partition S0 into 2 groups: those .lt. AM, those
C     .ge. AM IWRK(2*N+I)= how many pairs (X(I)+X(J)) in row I
C     less than AM
C
      N2 = 2*N
      J = N
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
C        Have we hit the diagonal ?
C
         IF (J.LT.I) THEN
            GO TO 40
C
C           Shall we move left ?
C
         ELSE IF (X(I)+X(J).GE.AM) THEN
            J = J - 1
            GO TO 20
         END IF
C
C           We're done in this row
C
         IWRK(N2+I) = J - I + 1
C
C           SQ = Total number of pairs less than AM
C
         SQ = SQ + IWRK(N2+I)
C
   40 CONTINUE
      RETURN
      END
