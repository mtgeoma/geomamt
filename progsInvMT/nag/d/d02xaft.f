      SUBROUTINE D02XAF(XSOL,X,COUT,N,Y,W,IW,SOL,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     CALCULATES SOLUTION BETWEEN STEPS OF D02PAF OR
C     D02KDY USING QUINTIC HERMITE INTERPOLATION
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02XAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, XSOL
      INTEGER           IFAIL, IW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  COUT(14), SOL(N), W(IW,5), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  H1, H2, S, S2, SP1, SP12, T, T2, TPS, TPS2,
     *                  TSP1, TSP12, TT, TT2
      INTEGER           I, IND
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Executable Statements ..
      IND = 0
      IF (IW.GE.N .AND. N.GT.0) GO TO 20
C     DIMENSIONS WRONG
      IND = 1
      GO TO 280
   20 H1 = X - COUT(4)
      H2 = COUT(4) - COUT(5)
      IF (H1.EQ.0.0D0 .OR. H2.EQ.0.0D0) GO TO 40
      IF (SIGN(1.D0,H1).EQ.SIGN(1.D0,H2)) GO TO 60
   40 CONTINUE
C     INTERPOLATION POINTS NOT ORDERED
      IND = 2
      GO TO 280
   60 IF (XSOL.EQ.X) GO TO 160
      IF (XSOL.EQ.COUT(4)) GO TO 200
      IF (XSOL.EQ.COUT(5)) GO TO 240
      IF (SIGN(1.D0,X-XSOL).EQ.SIGN(1.D0,XSOL-COUT(5))) GO TO 80
C     XSOL NOT AN INTERPOLATION POINT
      IND = 3
   80 S = H1/H2
      S2 = S*S
      SP1 = S + 1.D0
      SP12 = SP1*SP1
      T = (XSOL-COUT(4))
      IF (SIGN(1.D0,T).EQ.SIGN(1.D0,H1)) GO TO 120
      T = -T/H2
      T2 = T*T
      TT = 1.D0 - T
      TT2 = TT*TT
      TPS = T + S
      TPS2 = TPS*TPS
      DO 100 I = 1, N
         SOL(I) = T2*TPS2*(W(I,4)+TT*(H2*W(I,5)+2.D0*W(I,4)*(2.D0+S)
     *            /SP1))/SP12 + TT2*TPS2*(W(I,2)-T*(H2*W(I,3)
     *            -2.D0*W(I,2)*(S-1.D0)/S))/S2 + T2*TT2*(Y(I)
     *            -TPS*(H1*W(I,1)-2.D0*Y(I)*(2.D0*S+1.D0)/SP1)/S)
     *            /(S2*SP12)
  100 CONTINUE
      GO TO 280
  120 T = T/H1
      T2 = T*T
      TT = 1.D0 - T
      TT2 = TT*TT
      TSP1 = T*S + 1.D0
      TSP12 = TSP1*TSP1
      DO 140 I = 1, N
         SOL(I) = T2*TT2*S2*S2*(W(I,4)+TSP1*(H2*W(I,5)+2.D0*W(I,4)
     *            *(2.D0+S)/SP1))/SP12 + TSP12*TT2*(W(I,2)+T*(H1*W(I,3)
     *            -2.D0*W(I,2)*(S-1.D0))) + T2*TSP12*(Y(I)-TT*(H1*W(I,1)
     *            -2.D0*Y(I)*(2.D0*S+1.D0)/SP1))/SP12
  140 CONTINUE
      GO TO 280
  160 DO 180 I = 1, N
         SOL(I) = Y(I)
  180 CONTINUE
      GO TO 280
  200 DO 220 I = 1, N
         SOL(I) = W(I,2)
  220 CONTINUE
      GO TO 280
  240 DO 260 I = 1, N
         SOL(I) = W(I,4)
  260 CONTINUE
  280 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
