      SUBROUTINE G07EBU(N,X,M,Y,IWRK,AM,AMN,AMX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     This routine uses the midrange of set S as partition element
C     when ties are likely -- or to get the average of the last 2
C     elements
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AM, AMN, AMX
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(M)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      INTEGER           I, LBI, RBI
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      AMX = Y(1) - X(N)
      AMN = Y(M) - X(1)
      DO 20 I = 1, N
C
C        Skip this row if no element in it is in set S on this step.
C
         IF (IWRK(I).LE.IWRK(N+I)) THEN
            LBI = IWRK(I)
C
C           Get the smallest in this row.
C
            AMN = MIN(AMN,Y(LBI)-X(N-I+1))
            RBI = IWRK(N+I)
C
C           Get the largest in this row.
C
            AMX = MAX(AMX,Y(RBI)-X(N-I+1))
         END IF
   20 CONTINUE
C
      AM = (AMX+AMN)/2.D0
C
C     Be careful to cut off something -- roundoff can do wierd things
C
      IF (AM.LE.AMN .OR. AM.GT.AMX) AM = AMX
C
      RETURN
      END
