      SUBROUTINE D05BYN(IORDER,ALFA)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     -----------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++
C      This subroutine contains the BDF first
C      (reversed) charactrestic polynomial of orders
C      4 to 6.
C     +++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  ALFA(0:IORDER)
C     .. Local Scalars ..
      DOUBLE PRECISION  A12TH, A60TH
C     .. Executable Statements ..
      IF (IORDER.EQ.4) THEN
         A12TH = 1.D0/12.D0
         ALFA(0) = 25.D0*A12TH
         ALFA(1) = -48.D0*A12TH
         ALFA(2) = 36.D0*A12TH
         ALFA(3) = -16.D0*A12TH
         ALFA(4) = 3.D0*A12TH
C
      ELSE IF (IORDER.EQ.5) THEN
         A60TH = 1.D0/60.D0
         ALFA(0) = 137.D0*A60TH
         ALFA(1) = -300.D0*A60TH
         ALFA(2) = 300.D0*A60TH
         ALFA(3) = -200.D0*A60TH
         ALFA(4) = 75.D0*A60TH
         ALFA(5) = -12.D0*A60TH
C
      ELSE IF (IORDER.EQ.6) THEN
         A60TH = 1.D0/60.D0
         ALFA(0) = 147.D0*A60TH
         ALFA(1) = -360.D0*A60TH
         ALFA(2) = 450.D0*A60TH
         ALFA(3) = -400.D0*A60TH
         ALFA(4) = 225.D0*A60TH
         ALFA(5) = -72.D0*A60TH
         ALFA(6) = 10.D0*A60TH
C
      END IF
C
      RETURN
      END
