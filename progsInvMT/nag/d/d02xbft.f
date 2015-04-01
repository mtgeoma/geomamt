      SUBROUTINE D02XBF(XSOL,X,COUT,N,Y,W,IW,M,SOL,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     CALCULATES SOLUTION BETWEEN STEPS OF D02PAF OR
C     D02KDY USING QUINTIC HERMITE INTERPOLATION
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02XBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SOL, X, XSOL
      INTEGER           IFAIL, IW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  COUT(14), W(IW,5), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  H1, H2, S, S2, SP1, SP12, T, T2, TPS, TPS2,
     *                  TSP1, TSP12, TT, TT2
      INTEGER           IND
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Executable Statements ..
      IND = 0
      IF (IW.GE.N .AND. N.GT.0 .AND. M.GT.0 .AND. M.LE.N) GO TO 20
C     DIMENSIONS WRONG
      IND = 1
      GO TO 180
   20 H1 = X - COUT(4)
      H2 = COUT(4) - COUT(5)
      IF (H1.EQ.0.0D0 .OR. H2.EQ.0.0D0) GO TO 40
      IF (SIGN(1.D0,H1).EQ.SIGN(1.D0,H2)) GO TO 60
   40 CONTINUE
C     INTERPOLATION POINTS NOT ORDERED
      IND = 2
      GO TO 180
   60 IF (XSOL.EQ.X) GO TO 120
      IF (XSOL.EQ.COUT(4)) GO TO 140
      IF (XSOL.EQ.COUT(5)) GO TO 160
      IF (SIGN(1.D0,X-XSOL).EQ.SIGN(1.D0,XSOL-COUT(5))) GO TO 80
C     XSOL NOT AN INTERPOLATION POINT
      IND = 3
   80 S = H1/H2
      S2 = S*S
      SP1 = S + 1.D0
      SP12 = SP1*SP1
      T = (XSOL-COUT(4))
      IF (SIGN(1.D0,T).EQ.SIGN(1.D0,H1)) GO TO 100
      T = -T/H2
      T2 = T*T
      TT = 1.D0 - T
      TT2 = TT*TT
      TPS = T + S
      TPS2 = TPS*TPS
      SOL = T2*TPS2*(W(M,4)+TT*(H2*W(M,5)+2.D0*W(M,4)*(2.D0+S)/SP1))
     *      /SP12 + TT2*TPS2*(W(M,2)-T*(H2*W(M,3)-2.D0*W(M,2)*(S-1.D0)
     *      /S))/S2 + T2*TT2*(Y(M)-TPS*(H1*W(M,1)-2.D0*Y(M)
     *      *(2.D0*S+1.D0)/SP1)/S)/(S2*SP12)
      GO TO 180
  100 T = T/H1
      T2 = T*T
      TT = 1.D0 - T
      TT2 = TT*TT
      TSP1 = T*S + 1.D0
      TSP12 = TSP1*TSP1
      SOL = T2*TT2*S2*S2*(W(M,4)+TSP1*(H2*W(M,5)+2.D0*W(M,4)*(2.D0+S)
     *      /SP1))/SP12 + TSP12*TT2*(W(M,2)+T*(H1*W(M,3)-2.D0*W(M,2)
     *      *(S-1.D0))) + T2*TSP12*(Y(M)-TT*(H1*W(M,1)-2.D0*Y(M)
     *      *(2.D0*S+1.D0)/SP1))/SP12
      GO TO 180
  120 SOL = Y(M)
      GO TO 180
  140 SOL = W(M,2)
      GO TO 180
  160 SOL = W(M,4)
  180 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
