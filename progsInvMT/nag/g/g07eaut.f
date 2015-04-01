      SUBROUTINE G07EAU(N,X,IWRK,AM,AMN,AMX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Auxillary routine for G07EAW which uses the exact method for
C     computing the Hodges-Lehmann estimator and the corresponding
C     confidence interval.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AM, AMN, AMX
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      INTEGER           I, LBI, RBI
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Use the midrange of set S as partition element when ties are
C       likely -- or get the average of the last 2 elements
C
      AMX = X(1) + X(1)
      AMN = X(N) + X(N)
      DO 20 I = 1, N
C
C        Skip this row if no element in it is in set S on this step
C
         IF (IWRK(I).LE.IWRK(N+I)) THEN
            LBI = IWRK(I)
C
C           Get the smallest in this row
C
            AMN = MIN(AMN,X(LBI)+X(I))
            RBI = IWRK(N+I)
C
C           Get the largest in this row
C
            AMX = MAX(AMX,X(RBI)+X(I))
         END IF
   20 CONTINUE
      AM = (AMX+AMN)/2.D0
C
C     Be careful to cut off something -- roundoff can do wierd things
C
      IF (AM.LE.AMN .OR. AM.GT.AMX) AM = AMX
      RETURN
      END
