      SUBROUTINE D05BWJ(ALFA,WT,STWT,LDSW,IORDER,LENWT)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Scalar Arguments ..
      INTEGER           IORDER, LDSW, LENWT
C     .. Array Arguments ..
      DOUBLE PRECISION  ALFA(0:IORDER-1), STWT(LDSW,0:IORDER-1),
     *                  WT(0:LENWT)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACSUM, FACTOR, SUM
      INTEGER           I, IORDM1, J, K
C     .. Local Arrays ..
      DOUBLE PRECISION  V(6), VNJ(0:4,0:4)
C     .. Executable Statements ..
C
      IORDM1 = IORDER - 1
      FACTOR = 1.D0/ALFA(0)
      WT(0) = FACTOR
      DO 20 I = 1, IORDM1 - 1
         V(I) = 1.D0
   20 CONTINUE
      V(IORDM1) = 1.D0 - WT(0)
C
C     ...  Start evaluting the convolution weights ...
C
      DO 80 I = 2, LENWT + 1
         SUM = 0.D0
         DO 40 J = 1, IORDM1
            SUM = SUM + ALFA(J)*V(IORDER-J)
   40    CONTINUE
         DO 60 J = 1, IORDER - 2
            V(J) = V(J+1)
   60    CONTINUE
         V(IORDM1) = -FACTOR*SUM
         WT(I-1) = 1.D0 - V(IORDM1)
   80 CONTINUE
C
C     ... Start the evaluating the starting weights  ...
C
      DO 120 J = 0, IORDM1
         SUM = 0.D0
         DO 100 I = 0, IORDM1 - 1
            SUM = SUM + ALFA(I)*STWT(IORDM1-I,J)
  100    CONTINUE
         V(J+1) = SUM
  120 CONTINUE
      VNJ(0,0) = -V(1)
      DO 160 I = 1, IORDM1
         VNJ(0,I) = -V(I+1)
         DO 140 J = 0, IORDM1
            VNJ(I,J) = STWT(I,J) - V(J+1)
  140    CONTINUE
  160 CONTINUE
C
      DO 240 K = IORDER, (LENWT+IORDER)
         DO 220 J = 0, IORDM1
            SUM = 0.D0
            DO 180 I = 1, IORDM1
               SUM = SUM + ALFA(I)*VNJ(IORDER-I,J)
  180       CONTINUE
            FACSUM = -FACTOR*SUM
            STWT(K,J) = FACSUM + V(J+1)
            DO 200 I = 1, IORDM1 - 1
               VNJ(I,J) = VNJ(I+1,J)
  200       CONTINUE
            VNJ(IORDM1,J) = FACSUM
  220    CONTINUE
  240 CONTINUE
C
      RETURN
      END
