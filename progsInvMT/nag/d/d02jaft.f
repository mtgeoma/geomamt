      SUBROUTINE D02JAF(N,CF,BC,X0,X1,K1,KP,C,W,LW,IW,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     BC
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02JAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X0, X1
      INTEGER           IFAIL, K1, KP, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(K1), W(LW)
      INTEGER           IW(K1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF
      EXTERNAL          CF
C     .. Subroutine Arguments ..
      EXTERNAL          BC
C     .. Scalars in Common ..
      DOUBLE PRECISION  Y0, Y1
      INTEGER           M2
C     .. Local Scalars ..
      INTEGER           IA1, J, LA, N1
C     .. Local Arrays ..
      INTEGER           L(1), M(1)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02JAY, D02JAZ, D02TGY, D02TGZ
C     .. Common blocks ..
      COMMON            /AD02JA/Y0, Y1
      COMMON            /BD02JA/M2
C     .. Executable Statements ..
      N1 = KP + N
      J = 1
      IF (N.LE.0 .OR. K1.LT.N+1 .OR. N1.LT.K1) GO TO 20
      IA1 = 2*K1 + 2
      LA = N1*IA1
      J = 2
      IF (LW.LT.LA+K1*7) GO TO 20
      M(1) = N
      L(1) = N
      Y0 = X0
      Y1 = X1
      M2 = N
      CALL D02TGZ(1,M,L,X0,X1,K1,KP,C,K1,D02JAY,D02JAZ,D02TGY,D02TGY,CF,
     *            BC,W,N1,IA1,W(LA+1),K1,IW,J)
      IF (J.NE.0) GO TO 20
      IFAIL = 0
      GO TO 40
   20 IFAIL = P01ABF(IFAIL,J,SRNAME,0,P01REC)
   40 RETURN
      END
