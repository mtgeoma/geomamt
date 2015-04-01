      SUBROUTINE D05BAL(NBEG,NEND,NSTART,VG,VK,WKY,WT,L1WT,INC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     -------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine EVAluates the LAG-trems involving
C     convolution weights.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -------------------------------------------------------------
C
C
C     .. Scalar Arguments ..
      INTEGER           INC, L1WT, NBEG, NEND, NSTART
C     .. Array Arguments ..
      DOUBLE PRECISION  VG(0:NBEG-1), VK(0:INC*NEND), WKY(0:NEND),
     *                  WT(0:L1WT)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUMLAG
      INTEGER           NJ, NN, NSTAP1
C     .. Executable Statements ..
C
      NSTAP1 = NSTART + 1
      DO 100 NN = NBEG, NEND
         SUMLAG = 0.D0
         IF (NN.LE.L1WT+NSTAP1) THEN
            DO 20 NJ = NSTAP1, NBEG - 1
               SUMLAG = SUMLAG + WT(NN-NJ)*VK(INC*(NN-NJ))*VG(NJ)
   20       CONTINUE
         ELSE
            IF (NN.LE.L1WT+NBEG-1) THEN
               DO 40 NJ = NSTAP1, NN - L1WT - 1
                  SUMLAG = SUMLAG + VK(INC*(NN-NJ))*VG(NJ)
   40          CONTINUE
               DO 60 NJ = NN - L1WT, NBEG - 1
                  SUMLAG = SUMLAG + WT(NN-NJ)*VK(INC*(NN-NJ))*VG(NJ)
   60          CONTINUE
            ELSE
               DO 80 NJ = NSTAP1, NBEG - 1
                  SUMLAG = SUMLAG + VK(INC*(NN-NJ))*VG(NJ)
   80          CONTINUE
            END IF
         END IF
         WKY(NN) = SUMLAG + WKY(NN)
  100 CONTINUE
      RETURN
      END
