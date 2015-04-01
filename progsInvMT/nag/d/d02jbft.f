      SUBROUTINE D02JBF(N,CF,BC,X0,X1,K1,KP,C,IC,W,LW,IW,LIW,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     BC
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02JBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X0, X1
      INTEGER           IC, IFAIL, K1, KP, LIW, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(IC,N), W(LW)
      INTEGER           IW(LIW)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF
      EXTERNAL          CF
C     .. Subroutine Arguments ..
      EXTERNAL          BC
C     .. Scalars in Common ..
      DOUBLE PRECISION  Y0, Y1
      INTEGER           NN
C     .. Local Scalars ..
      INTEGER           I, IA1, IRB, J, LA, N1, N2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02JBY, D02JBZ, D02TGY, D02TGZ
C     .. Common blocks ..
      COMMON            /AD02JB/Y0, Y1, NN
C     .. Executable Statements ..
      J = 1
      IF (N.LE.0 .OR. K1.LT.2 .OR. KP+1.LT.K1) GO TO 60
      N1 = N*KP + N
      IRB = N*K1
      IA1 = 2*IRB + 2
      LA = N1*IA1
      J = 2
      IF (LW.LT.LA+IRB*7 .OR. LIW.LT.IRB+2*N) GO TO 60
      N2 = N + N
      DO 20 I = 1, N2
         IW(I) = 1
   20 CONTINUE
      Y0 = X0
      Y1 = X1
      NN = N
      CALL D02TGZ(N,IW,IW(N+1),X0,X1,K1,KP,C,IC,D02JBY,D02JBZ,D02TGY,
     *            D02TGY,CF,BC,W,N1,IA1,W(LA+1),IRB,IW(2*N+1),J)
      IF (J.EQ.0) GO TO 40
      IFAIL = 0
      GO TO 60
   40 IFAIL = P01ABF(IFAIL,J,SRNAME,0,P01REC)
   60 RETURN
      END
