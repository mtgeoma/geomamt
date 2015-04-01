      SUBROUTINE D05ABF(K,G,LAMBDA,A,B,ODOREV,EV,N,CM,F1,WK,NMAX,NT2P1,
     *                  F,C,IFAIL)
C     MARK 6 RELEASE. NAG COPYRIGHT 1977.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NPL DNAC LIBRARY SUBROUTINE F2925.
C
C     THIS SUBROUTINE SOLVES THE NON-SINGULAR FREDHOLM INTEGRAL
C     EQUATION OF THE SECOND KIND-
C     F(X) - LAMBDA*(INTEGRAL FROM A TO B OF K(X,S)*F(S) DS) = G(X)
C     FOR A .LE. X .LE. B, BY THE METHOD OF S.E. EL-GENDI
C     (COMPUTER J., 1969, VOL.12, PP.282-287).
C
C     THE SUBROUTINE USES AUXILIARY SUBROUTINES F04AAF
C     AND FUNCTIONS X01AAF,P01AAF AND C06DBF.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D05ABF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, LAMBDA
      INTEGER           IFAIL, N, NMAX, NT2P1
      LOGICAL           EV, ODOREV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), CM(NMAX,NMAX), F(N), F1(NMAX,1),
     *                  WK(2,NT2P1)
C     .. Function Arguments ..
      DOUBLE PRECISION  G, K
      EXTERNAL          G, K
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, B1, C1, D1, D11, FT, NR, P, P2, PI, S, WJ,
     *                  WJ1, X
      INTEGER           I, IERROR, J, N1, N1MJP2, N1P1, N2, N2M1, N2P1,
     *                  N3, NN
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  C06DBF, X01AAF
      INTEGER           P01ABF
      EXTERNAL          C06DBF, X01AAF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE
C     .. Executable Statements ..
      A1 = 0.5D0*(A+B)
      B1 = 0.5D0*(B-A)
      C1 = B1*LAMBDA
      IF (B1.GT.0.0D0) GO TO 20
      IERROR = 1
      GO TO 460
   20 IERROR = 0
      NN = N - 1
      IF ( .NOT. ODOREV) GO TO 40
      IF ( .NOT. EV) GO TO 60
      N1 = 2*N - 1
      GO TO 80
   40 N1 = NN
      GO TO 80
   60 N1 = 2*N
   80 IF (N1.NE.0) GO TO 120
C
C     TRIVIAL CASE
C
      B1 = 1.0D0 - C1*K(A1,A1)
      IF (B1.NE.0.0D0) GO TO 100
      IERROR = 2
      GO TO 460
  100 F(1) = G(A1)/B1
      C(1) = 2.0D0*F(1)
      GO TO 480
C
C     GENERAL CASE
C
  120 PI = X01AAF(PI)
      N2 = N1/2
      N3 = 2*N2
      NR = 1.0D0/DBLE(N1)
      P = PI*NR
C
C     FORM WEIGHTS FOR QUADRATURE
C
      FT = 1.0D0/DBLE(N3*N3-1)
      WK(1,1) = DBLE(2*N3-N1)*NR*FT
      P2 = 2.0D0*P
      IF (N2.EQ.0) GO TO 200
      N2M1 = N2 - 1
      DO 180 J = 1, N2
         WJ = 0.5D0
         IF (N2M1.EQ.0) GO TO 160
         DO 140 I = 1, N2M1
            WJ = WJ + COS(DBLE(I*J)*P2)/DBLE(1-4*I*I)
  140    CONTINUE
  160    WJ1 = FT*COS(DBLE(N2*J)*P2)
         IF (N3.EQ.N1) WJ = WJ - 0.5D0*WJ1
         IF (N3.NE.N1) WJ = WJ - WJ1
         WK(1,J+1) = 4.0D0*WJ*NR
  180 CONTINUE
  200 N2P1 = N2 + 1
      IF (ODOREV) GO TO 240
      DO 220 J = 1, N2P1
         N1MJP2 = N1 - J + 2
         WK(1,N1MJP2) = WK(1,J)
  220 CONTINUE
C
C     FORM AND SOLVE EQUATIONS
C
  240 N1P1 = N1 + 1
      DO 260 I = 1, N1P1
         WK(2,I) = A1 + B1*COS(DBLE(I-1)*P)
  260 CONTINUE
      DO 320 I = 1, N
         X = WK(2,I)
         F1(I,1) = G(X)
         DO 300 J = 1, N
            S = WK(2,J)
            D1 = WK(1,J)*K(X,S)
            IF ( .NOT. ODOREV .OR. J.GT.N2P1) GO TO 280
            N1MJP2 = N1 - J + 2
            S = WK(2,N1MJP2)
            D11 = WK(1,J)*K(X,S)
            IF (EV) D1 = D1 + D11
            IF ( .NOT. EV) D1 = D1 - D11
  280       CM(I,J) = -C1*D1
            IF (I.EQ.J) CM(I,J) = 1.0D0 + CM(I,J)
  300    CONTINUE
  320 CONTINUE
      I = 1
      CALL F04AAF(CM,NMAX,F1,NMAX,N,1,F1,NMAX,F,I)
      IF (I.EQ.0) GO TO 340
      IERROR = 2
      GO TO 460
  340 DO 360 I = 1, N
         F(I) = F1(I,1)
  360 CONTINUE
C
C     FIND CHEBYSHEV COEFFICIENTS OF SOLUTION
C
      IF ( .NOT. ODOREV) F(N) = F(N)*0.5D0
      DO 440 I = 1, N
         IF ( .NOT. ODOREV) GO TO 400
         IF ( .NOT. EV) GO TO 380
         X = COS(2.0D0*DBLE(I-1)*P)
         GO TO 420
  380    X = COS((2.0D0*DBLE(I-1)+1.0D0)*P)
         GO TO 420
  400    X = COS(DBLE(I-1)*P)
  420    C(I) = C06DBF(X,F,N,1)*NR*2.0D0
         IF (ODOREV) C(I) = 2.0D0*C(I)
  440 CONTINUE
      IF (ODOREV) GO TO 480
      C(N) = 0.5D0*C(N)
      F(N) = 2.0D0*F(N)
      GO TO 480
  460 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      GO TO 500
  480 IFAIL = 0
  500 RETURN
      END
