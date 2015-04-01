      SUBROUTINE D05BWN(IORDER,LENWT,ALFA,WT,STWT,LDSW)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     -----------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++
C      This subroutine generates the BDF weights of
C      orders 2 to 5.
C     +++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IORDER, LDSW, LENWT
C     .. Array Arguments ..
      DOUBLE PRECISION  ALFA(0:IORDER-1), STWT(LDSW,0:IORDER-1),
     *                  WT(0:LENWT)
C     .. Local Scalars ..
      DOUBLE PRECISION  A12TH, A60TH, A720TH, ATHIRD, TRDTI
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          D05BWJ
C     .. Executable Statements ..
C
      IF (IORDER.EQ.2) THEN
C
C        ----   Start the generation of  weights associated  ----
C        ----   with the BDF2 formula                        ----
C
         ATHIRD = 1.D0/3.D0
         WT(0) = 2.D0*ATHIRD
         STWT(1,0) = 0.5D0
         STWT(1,1) = 0.5D0
C
         DO 40 I = 2, LENWT + 1
            TRDTI = 1.D0 - ATHIRD**I
            WT(I-1) = TRDTI
C
            DO 20 J = 0, 1
               STWT(I,J) = 0.75D0*TRDTI
   20       CONTINUE
C
   40    CONTINUE
C
         DO 80 I = LENWT + 2, LENWT + IORDER
            TRDTI = 1.D0 - ATHIRD**I
            DO 60 J = 0, 1
               STWT(I,J) = 0.75D0*TRDTI
   60       CONTINUE
   80    CONTINUE
C
C
      ELSE
C
C        ----   End of BDF2   ----
C        ----   Start the generation of  weights associated  ----
C        ----   with the BDF3 formula                        ----
C
         IF (IORDER.EQ.3) THEN
            A12TH = 1.D0/12.D0
            ALFA(0) = 11.D0/6.D0
            ALFA(1) = -7.D0/6.D0
            ALFA(2) = 1.D0/3.D0
            STWT(1,0) = 5.D0*A12TH
            STWT(1,1) = 8.D0*A12TH
            STWT(1,2) = -A12TH
            STWT(2,0) = 4.D0*A12TH
            STWT(2,1) = 16.D0*A12TH
            STWT(2,2) = 4.D0*A12TH
C
C           ----   End BDF3  ----
C           ----   Start the generation of  weights associated  ----
C           ----   with the BDF4 formula                        ----
C
         ELSE IF (IORDER.EQ.4) THEN
            A12TH = 1.D0/12.D0
            ALFA(0) = 25.D0*A12TH
            ALFA(1) = -23.D0*A12TH
            ALFA(2) = 13.D0*A12TH
            ALFA(3) = -.25D0
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
C
C           ----   End  BDF4  ----
C           ----   Start the generation of  weights associated  ----
C           ----   with the BDF5 formula                        ----
C
         ELSE IF (IORDER.EQ.5) THEN
            A60TH = 1.D0/60.D0
            A720TH = 1.D0/720.D0
            ALFA(0) = 137.D0*A60TH
            ALFA(1) = -163.D0*A60TH
            ALFA(2) = 137.D0*A60TH
            ALFA(3) = -63.D0*A60TH
            ALFA(4) = 0.2D0
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
         END IF
C
C        ----  End  BDF5  ----
C
         CALL D05BWJ(ALFA,WT,STWT,LDSW,IORDER,LENWT)
C
      END IF
      RETURN
      END
