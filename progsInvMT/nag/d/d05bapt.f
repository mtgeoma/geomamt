      SUBROUTINE D05BAP(CG,CF,ALIM,WVK,WVG,WKY,VF,VG,H,THETA,B,A,IQ,
     *                  NSTART)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine computes the solution of a convolution
C     Volterra equation, using an explicit Runge-Kutta scheme of SDP.
C     The purpose of this routine, merely, is to compute the starting
C     values.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, H
      INTEGER           IQ, NSTART
C     .. Array Arguments ..
      DOUBLE PRECISION  A(6,0:5), B(0:IQ-1), THETA(0:IQ-1),
     *                  VF(0:NSTART), VG(0:NSTART), WKY(0:NSTART),
     *                  WVG(0:NSTART,0:IQ-1), WVK(0:NSTART-1,0:IQ-1,
     *                  0:IQ-1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF, CG
      EXTERNAL          CF, CG
C     .. Local Scalars ..
      DOUBLE PRECISION  SUMAUX, SUMLAG, TN, TNR
      INTEGER           I, IQM1, IS, J, NJ, NN, NN1
C     .. Executable Statements ..
C
      IQM1 = IQ - 1
      TN = ALIM
      WVG(0,0) = CG(TN,VF(0))
      DO 100 NN = 0, NSTART - 1
         DO 80 NJ = 1, IQM1
            TNR = TN + H*THETA(NJ)
C
C           ---- Evaluate the lag-term. ----
C
            SUMLAG = 0.D0
            DO 40 IS = 0, IQM1
               SUMAUX = 0.D0
               DO 20 I = 0, NN - 1
                  SUMAUX = SUMAUX + WVK(NN-I,NJ,IS)*WVG(I,IS)
   20          CONTINUE
               SUMLAG = SUMLAG + B(IS)*SUMAUX
   40       CONTINUE
            DO 60 J = 0, NJ - 1
               SUMLAG = SUMLAG + A(NJ,J)*WVK(0,NJ,J)*WVG(NN,J)
   60       CONTINUE
            SUMLAG = H*SUMLAG + CF(TNR)
            WVG(NN,NJ) = CG(TNR,SUMLAG)
   80    CONTINUE
         NN1 = NN + 1
         TN = TN + H
         WKY(NN1) = SUMLAG
         WVG(NN1,0) = WVG(NN,IQM1)
         VG(NN1) = H*WVG(NN1,0)
  100 CONTINUE
      RETURN
      END
