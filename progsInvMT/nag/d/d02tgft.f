      SUBROUTINE D02TGF(N,M,L,X0,X1,K1,KP,C,IC,COEFF,BDYC,W,LW,IW,LIW,
     *                  IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     BDYC,COEFF
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02TGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X0, X1
      INTEGER           IC, IFAIL, K1, KP, LIW, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(IC,N), W(LW)
      INTEGER           IW(LIW), L(N), M(N)
C     .. Subroutine Arguments ..
      EXTERNAL          BDYC, COEFF
C     .. Local Scalars ..
      INTEGER           I, IA1, IER, IRB, LA, N1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02TGX
      INTEGER           P01ABF
      EXTERNAL          D02TGX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02TGV, D02TGW, D02TGY, D02TGZ
C     .. Executable Statements ..
      IER = 1
      IF (N.LE.0) GO TO 40
      N1 = N*KP
      DO 20 I = 1, N
         IF (L(I).LT.0) GO TO 40
         N1 = N1 + L(I)
   20 CONTINUE
      IRB = N*K1
      IA1 = 2*IRB + 2
      LA = N1*IA1
      IER = 2
      IF (LW.LT.LA+IRB*7 .OR. LIW.LT.IRB) GO TO 40
      CALL D02TGZ(N,M,L,X0,X1,K1,KP,C,IC,D02TGV,D02TGW,COEFF,BDYC,
     *            D02TGX,D02TGY,W,N1,IA1,W(LA+1),IRB,IW,IER)
      IF (IER.NE.0) GO TO 40
      IFAIL = 0
      GO TO 60
   40 IFAIL = P01ABF(IFAIL,IER,SRNAME,0,P01REC)
   60 RETURN
      END
