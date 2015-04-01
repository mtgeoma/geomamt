      SUBROUTINE D05BWM(IORDER,WT,STWT)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ---------------------------------------
C     +++++++++++++++++++++++++++++++++++++++
C      This subroutine generates the weights
C      associated with Adam's formulae of
C      the orders 3 to 6.
C     +++++++++++++++++++++++++++++++++++++++
C     ---------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  STWT(9,0:4), WT(0:5)
C     .. Local Scalars ..
      DOUBLE PRECISION  A12TH, A720TH
      INTEGER           I
C     .. Executable Statements ..
C
      IF (IORDER.EQ.3) THEN
         A12TH = 1.D0/12.D0
         WT(0) = 5.D0*A12TH
         WT(1) = 13.D0*A12TH
         WT(2) = 1.D0
         STWT(1,0) = 0.5D0
         STWT(1,1) = 0.5D0
         STWT(2,0) = 5.D0*A12TH
         STWT(2,1) = 14.D0*A12TH
         STWT(3,0) = WT(0)
         STWT(3,1) = WT(1)
C
      ELSE IF (IORDER.EQ.4) THEN
         A12TH = 1.D0/12.D0
         WT(0) = 4.5D0*A12TH
         WT(1) = 14.D0*A12TH
         WT(2) = 11.5D0*A12TH
         WT(3) = 1.D0
         STWT(1,0) = 5.D0*A12TH
         STWT(1,1) = 8.D0*A12TH
         STWT(1,2) = -A12TH
         STWT(2,0) = 4.D0*A12TH
         STWT(2,1) = 16.D0*A12TH
         STWT(2,2) = 4.D0*A12TH
         STWT(3,0) = 4.5D0*A12TH
         STWT(3,1) = 13.5D0*A12TH
         STWT(3,2) = 13.5D0*A12TH
         STWT(4,0) = 4.5D0*A12TH
         STWT(4,1) = 14.D0*A12TH
         STWT(4,2) = 11.D0*A12TH
         STWT(5,0) = WT(0)
         STWT(5,1) = WT(1)
         STWT(5,2) = WT(2)
C
      ELSE IF (IORDER.EQ.5) THEN
         A720TH = 1.D0/720.D0
         A12TH = 1.D0/12.D0
         WT(0) = 251.D0*A720TH
         WT(1) = 897.D0*A720TH
         WT(2) = 633.D0*A720TH
         WT(3) = 739.D0*A720TH
         WT(4) = 1.D0
         STWT(1,0) = 4.5D0*A12TH
         STWT(1,1) = 9.5D0*A12TH
         STWT(1,2) = -2.5D0*A12TH
         STWT(1,3) = 0.5D0*A12TH
         STWT(2,0) = 4.D0*A12TH
         STWT(2,1) = 16.D0*A12TH
         STWT(2,2) = 4.D0*A12TH
         STWT(2,3) = 0.D0
         STWT(3,0) = 4.5D0*A12TH
         STWT(3,1) = 13.5D0*A12TH
         STWT(3,2) = 13.5D0*A12TH
         STWT(3,3) = 4.5D0*A12TH
         STWT(4,0) = 251.D0*A720TH
         STWT(4,1) = 916.D0*A720TH
         STWT(4,2) = 546.D0*A720TH
         STWT(4,3) = 916.D0*A720TH
         STWT(5,0) = 251.D0*A720TH
         STWT(5,1) = 897.D0*A720TH
         STWT(5,2) = 652.D0*A720TH
         STWT(5,3) = 652.D0*A720TH
         STWT(6,0) = 251.D0*A720TH
         STWT(6,1) = 897.D0*A720TH
         STWT(6,2) = 633.D0*A720TH
         STWT(6,3) = 758.D0*A720TH
         STWT(7,0) = WT(0)
         STWT(7,1) = WT(1)
         STWT(7,2) = WT(2)
         STWT(7,3) = WT(3)
C
      ELSE IF (IORDER.EQ.6) THEN
         A720TH = 1.D0/720.D0
         WT(0) = 237.5D0*A720TH
         WT(1) = 951.D0*A720TH
         WT(2) = 552.D0*A720TH
         WT(3) = 793.D0*A720TH
         WT(4) = 706.5D0*A720TH
         WT(5) = 1.D0
         STWT(1,0) = 251.D0*A720TH
         STWT(1,1) = 646.D0*A720TH
         STWT(1,2) = -264.D0*A720TH
         STWT(1,3) = 106.D0*A720TH
         STWT(1,4) = -19.D0*A720TH
         STWT(2,0) = 232.D0*A720TH
         STWT(2,1) = 992.D0*A720TH
         STWT(2,2) = 192.D0*A720TH
         STWT(2,3) = 32.D0*A720TH
         STWT(2,4) = -8.D0*A720TH
         STWT(3,0) = 243.D0*A720TH
         STWT(3,1) = 918.D0*A720TH
         STWT(3,2) = 648.D0*A720TH
         STWT(3,3) = 378.D0*A720TH
         STWT(3,4) = -27.D0*A720TH
         STWT(4,0) = 224.D0*A720TH
         STWT(4,1) = 1024.D0*A720TH
         STWT(4,2) = 384.D0*A720TH
         STWT(4,3) = 1024.D0*A720TH
         STWT(4,4) = 224.D0*A720TH
         STWT(5,0) = 237.5D0*A720TH
         STWT(5,1) = 937.5D0*A720TH
         STWT(5,2) = 625.D0*A720TH
         STWT(5,3) = 625.D0*A720TH
         STWT(5,4) = 937.5D0*A720TH
         STWT(6,0) = 237.5D0*A720TH
         STWT(6,1) = 951.D0*A720TH
         STWT(6,2) = 538.5D0*A720TH
         STWT(6,3) = 866.D0*A720TH
         STWT(6,4) = 538.5D0*A720TH
         DO 20 I = 0, 2
            STWT(7,I) = WT(I)
   20    CONTINUE
         STWT(7,3) = 779.5D0*A720TH
         STWT(7,4) = STWT(7,3)
         DO 40 I = 0, 3
            STWT(8,I) = WT(I)
   40    CONTINUE
         STWT(8,4) = 693.D0*A720TH
         DO 60 I = 0, 4
            STWT(9,I) = WT(I)
   60    CONTINUE
      END IF
C
      RETURN
      END
