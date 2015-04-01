      SUBROUTINE D05BAS(CG,H,ALIM,WT,L1WT,VK,VG,WKY,INIT,NE,NSTART,INC,
     *                  INFO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     --------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine computes the SOLution of the convolution
C     VOLterra equation. WKY(i) approximates the true solution
C     y(t), t = i*H, i = INIT + 1, . . . , NE. At each stage
C     the value of  H*G(i*H, WKY(i)) is stored in VG(i).
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     --------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, H
      INTEGER           INC, INFO, INIT, L1WT, NE, NSTART
C     .. Array Arguments ..
      DOUBLE PRECISION  VG(0:NE), VK(0:INC*NE), WKY(0:NE), WT(0:L1WT)
C     .. Function Arguments ..
      DOUBLE PRECISION  CG
      EXTERNAL          CG
C     .. Scalars in Common ..
      DOUBLE PRECISION  GTOL
C     .. Local Scalars ..
      DOUBLE PRECISION  BL, BU, GEVAL, GFUN, HS, SUMLAG, TN, WK0, Y0, Z0
      INTEGER           IFLAG, IND, INITP1, IR, NJ, NN
C     .. Local Arrays ..
      DOUBLE PRECISION  C(17)
C     .. External Subroutines ..
      EXTERNAL          C05AVF, C05AZF
C     .. Common blocks ..
      COMMON            /AD05BA/GTOL
C     .. Executable Statements ..
C
      INITP1 = INIT + 1
      WK0 = H*WT(0)*VK(0)
      TN = ALIM + INIT*H
      DO 160 NN = INITP1, NE
         TN = TN + H
C
C        ----  Start the evaluation of the lag-term.  ----
C
         SUMLAG = 0.D0
         IF (NN.LE.L1WT+INITP1) THEN
            DO 20 NJ = INITP1, NN - 1
               SUMLAG = SUMLAG + WT(NN-NJ)*VK(INC*(NN-NJ))*VG(NJ)
   20       CONTINUE
         ELSE
            DO 40 NJ = INITP1, NN - L1WT - 1
               SUMLAG = SUMLAG + VK(INC*(NN-NJ))*VG(NJ)
   40       CONTINUE
            DO 60 NJ = NN - L1WT, NN - 1
               SUMLAG = SUMLAG + WT(NN-NJ)*VK(INC*(NN-NJ))*VG(NJ)
   60       CONTINUE
         END IF
         SUMLAG = SUMLAG + WKY(NN)
C
C        ----  Solve the nonlinear equation.  ----
C
         IND = 1
         IFLAG = 1
         Y0 = WKY(NN-1)
         HS = 0.1D0
         BL = Y0 - 256.0D0*HS
         BU = Y0 + 256.0D0*HS
   80    CALL C05AVF(Y0,GFUN,HS,BL,BU,Z0,C,IND,IFLAG)
         IF (IND.EQ.0) GO TO 100
         GEVAL = CG(TN,Y0)
         GFUN = Y0 - WK0*GEVAL - SUMLAG
         GO TO 80
  100    CONTINUE
         IF (IFLAG.GT.0) THEN
            INFO = 4
            RETURN
         END IF
         IR = 1
         IFLAG = 1
         IND = 1
  120    CALL C05AZF(Y0,Z0,GFUN,GTOL,IR,C,IND,IFLAG)
         GEVAL = CG(TN,Y0)
         IF (IND.EQ.0) GO TO 140
         GFUN = Y0 - WK0*GEVAL - SUMLAG
         GO TO 120
  140    CONTINUE
         IF (IFLAG.GT.0) THEN
            INFO = 4
            RETURN
         END IF
         VG(NN) = H*GEVAL
         WKY(NN) = Y0
  160 CONTINUE
      RETURN
      END
