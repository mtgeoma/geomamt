      SUBROUTINE D02AGZ(X,Y,G,T,N,N1,M,H0,H,AUX,WSPACE,WSPAC1,PARAM)
C
C     USES THE LOGIC OF NAG LIBRARY ROUTINE D02ABF
C     EPS AND DUM MAY BE USED DOUBLE LENGTH.  EPS IS THE
C     SMALLEST REAL SUCH THAT 1+EPS>EPS.  SMAX
C     IS THE LARGEST INTEGER. DUM,ERR AND HS MAYBE DECLARED
C     DOUBLE PRECISION
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, H0, X
      INTEGER           M, N, N1, T
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), PARAM(N1), WSPAC1(N), WSPACE(N,9), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          AUX
C     .. Local Scalars ..
      DOUBLE PRECISION  DUM, EPS, ERR, HS, P, Q
      INTEGER           D, I, J, S, SMAX
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           X02BBF
      EXTERNAL          X02AJF, X02BBF
C     .. External Subroutines ..
      EXTERNAL          D02AGY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, AINT, MAX, DBLE, INT
C     .. Executable Statements ..
      EPS = X02AJF()
      SMAX = X02BBF(EPS)
      M = 0
      DUM = MAX(ABS(X),ABS(X+H0))
      IF (ABS(H0).LE.EPS*DUM) RETURN
      DUM = H/H0
      IF (ABS(DUM).GT.1.0D0 .OR. H.EQ.0.0D0) GO TO 20
      DUM = ABS(H0/H+0.9D0)
      IF (DUM.GT.DBLE(SMAX)) GO TO 160
      H = H0/AINT(DUM)
      GO TO 40
   20 H = H0
   40 P = 1.0D0
      IF (T.EQ.3) P = 0.0D0
      Q = 1.0D0
      IF (T.EQ.2) Q = 0.0D0
      DUM = H0/H + 0.1D0
      S = INT(DUM)
      HS = 1.0D-4*ABS(H)
   60 DO 80 I = 1, N
         WSPACE(I,9) = Y(I)
   80 CONTINUE
      CALL D02AGY(Y,X,H,N,N1,AUX,WSPACE,WSPAC1,PARAM)
      D = 0
      DO 100 I = 1, N
         ERR = G(I)*(P+Q*ABS(Y(I)))
         IF (WSPACE(I,2).GT.ERR) GO TO 120
         IF (40.0D0*WSPACE(I,2).GT.ERR) D = 1
  100 CONTINUE
      S = S - 1
      IF (S.EQ.0) RETURN
      DUM = DBLE(S)/2.0D0 + 0.1D0
      J = INT(DUM)*2
      IF (D.NE.0 .OR. J.NE.S) GO TO 60
      H = 2.0D0*H
      S = S/2
      GO TO 60
  120 X = X - H
      DO 140 I = 1, N
         Y(I) = WSPACE(I,9)
  140 CONTINUE
      IF (S.GT.SMAX/2) GO TO 160
      S = 2*S
      H = 0.5D0*H
      IF (ABS(H).GT.HS) GO TO 60
  160 M = 1
      RETURN
      END
