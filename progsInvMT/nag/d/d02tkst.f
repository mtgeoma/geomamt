      INTEGER FUNCTION D02TKS(A,XI,N)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C******************************************************************
C
C   Purpose:
C      to determine the interval in XI where A is located
C
C   Arguments:
C      A      - the point of interest
C      XI     - the current mesh  XI(i) < XI(i+1)
C      N      - the number of intervals in the mesh
C               (ie. the number of points in XI - 1)
C
C   Author:
C      R.W. Brankin, NAG Ltd, Aug 1994
C      (replaces linear search used in old COLNEW routine APPROX)
C
C******************************************************************
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        A
      INTEGER                 N
C     .. Array Arguments ..
      DOUBLE PRECISION        XI(N+1)
C     .. Local Scalars ..
      INTEGER                 I, IL, IR, J
C     .. Executable Statements ..
C
      J = 0
C
C Check for end points
C
      IF (A.EQ.XI(1)) THEN
         IL = 1
      ELSE IF (A.EQ.XI(N+1)) THEN
         IL = N
      ELSE IF (A.LT.XI(1) .OR. A.GT.XI(N+1)) THEN
         IL = -1
      ELSE
C
C Narrow down interval by bisection on interval number
C
         IL = 1
         IR = N + 1
   20    CONTINUE
         J = J + 1
         I = (IR+IL)/2
         IF (A.LT.XI(I)) THEN
            IR = I
         ELSE
            IL = I
         END IF
         IF (IR-IL.NE.1) GO TO 20
      END IF
      D02TKS = IL
      RETURN
      END
