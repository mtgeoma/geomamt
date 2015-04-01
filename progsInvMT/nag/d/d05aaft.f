      SUBROUTINE D05AAF(LAMBDA,A,B,K1,K2,G,F,C,N,IND,W1,W2,WD,NMAX,MN,
     *                  IFAIL)
C     MARK 5 RELEASE  NAG COPYRIGHT 1976
C     MARK 6 REVISED
C     MARK 9 REVISED. IER-304 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     BASED UPON NPL DNAC LIBRARY SUBROUTINE F2926
C     THIS SUBROUTINE SOLVES THE FREDHOLM INTEGRAL EQUATION OF THE
C     SECOND KIND
C     F(X) - LAMBDA * (INTEGRAL FROM A TO B OF K(X,S) * F(S) DS) =
C     G(X)
C     FOR A.LE.X.LE.B, WHEN THE KERNEL IS DEFINED IN TWO PARTS
C     K = K1 FOR A.LE.S.LE.X AND K = K2 FOR X.LT.S.LE.B.
C     THE METHOD USED IS THAT OF EL-GENDI, WHICH REQUIRES
C     THAT EACH OF THE FUNCTIONS K1 AND K2 SHOULD BE SMOOTH
C     AND NON SINGULAR FOR ALL X AND S IN THE CLOSED INTERVAL
C     (A,B).
C     THE SUBROUTINE USES AUXILIARY SUBROUTINE F04AAF
C     AND FUNCTIONS P01AAF AND C06DBF.
C     NAG COPYRIGHT 1976
C     MARK 5 RELEASE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D05AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, LAMBDA
      INTEGER           IFAIL, IND, MN, N, NMAX
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), F(N), W1(NMAX,MN), W2(MN,4), WD(MN)
C     .. Function Arguments ..
      DOUBLE PRECISION  G, K1, K2
      EXTERNAL          G, K1, K2
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, B1, C1, NR, P, PI, S, SIG, T0, T1, T2, V0,
     *                  V1, V2, X, X2
      INTEGER           I, I1, I2, ICM2, IDKS, IERROR, IGI1, IGI2, J,
     *                  J2, JDKS, M, M1, M2, MMIP1, MMJP1, N1, N1MIP2,
     *                  N1MJP2, N2, N2P1, N2T2
      LOGICAL           GEV, GEVOD, JODD, KCS
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
      KCS = .TRUE.
      GEVOD = .TRUE.
      GEV = .TRUE.
      IF ((IND.LE.0) .OR. (IND.GT.3)) KCS = .FALSE.
      IF (IND.EQ.1) GEV = .FALSE.
      IF (IND.EQ.3) GEVOD = .FALSE.
      ICM2 = NMAX + 1
      IGI1 = ICM2 + NMAX
      IGI2 = IGI1 + 1
      A1 = 0.5D0*(A+B)
      B1 = 0.5D0*(B-A)
      C1 = B1*LAMBDA
      IF (B1.GT.0.0D0) GO TO 20
      IERROR = 1
      GO TO 1260
   20 IERROR = 0
      IF (KCS .AND. GEVOD) GO TO 40
      N1 = N - 1
      GO TO 80
   40 IF (GEV) GO TO 60
      N1 = 2*N
      GO TO 80
   60 N1 = 2*N - 1
   80 IF (N1.NE.0) GO TO 120
C     TRIVIAL CASE
      P = 1.0D0 - C1*K1(A1,A1)
      IF (P.NE.0.0D0) GO TO 100
      GO TO 1240
  100 F(1) = G(A1)/P
      C(1) = 2.0D0*F(1)
C     END TRIVIAL CASE
      GO TO 1260
C     GENERAL CASE
  120 PI = X01AAF(PI)
      N2 = N1/2
      N2P1 = N2 + 1
      NR = 1.0D0/DBLE(N1)
      P = PI*NR
      M = N1 + 1
      IF (KCS) GO TO 140
      M1 = M
      GO TO 180
  140 IF ( .NOT. GEVOD .OR. GEV) GO TO 160
      M1 = 1
      GO TO 180
  160 M1 = N2 + 1
  180 IF (KCS .AND. .NOT. (GEVOD .AND. GEV)) GO TO 200
      M2 = 1
      GO TO 220
  200 M2 = N1 - N2
C     CHEBYSHEV POINTS AND RHS
  220 DO 240 I = 1, M
         W2(I,2) = COS(DBLE(I-1)*P)
         X = A1 + B1*W2(I,2)
         W2(I,4) = G(X)
  240 CONTINUE
      IF ( .NOT. KCS) GO TO 320
      IF (GEVOD .AND. .NOT. GEV) GO TO 280
      DO 260 I = 1, M1
         MMIP1 = M - I + 1
         W1(I,IGI1) = W2(I,4) + W2(MMIP1,4)
  260 CONTINUE
  280 IF (GEVOD .AND. GEV) GO TO 360
      DO 300 I = 1, M2
         MMIP1 = M - I + 1
         W1(I,IGI2) = W2(I,4) - W2(MMIP1,4)
  300 CONTINUE
      GO TO 360
  320 DO 340 I = 1, M1
         W1(I,IGI1) = W2(I,4)
  340 CONTINUE
C     FIRST ROW OF WEIGHTS
  360 DO 380 I = 1, N2P1
         WD(I) = 1.0D0/DBLE(1-4*(I-1)*(I-1))
  380 CONTINUE
      N2T2 = N2*2
      IF (N2T2.NE.N1) GO TO 400
      WD(N2P1) = 0.5D0*WD(N2P1)
  400 DO 420 J = 1, N2P1
         W2(J,3) = 4.0D0*NR*C06DBF(W2(J,2),WD,N2P1,2)
         N1MJP2 = N1 - J + 2
         W2(N1MJP2,3) = W2(J,3)
  420 CONTINUE
      W2(1,3) = 0.5D0*W2(1,3)
      W2(M,3) = W2(1,3)
C     FIRST ROW OF MATRIX
      X = B
      DO 440 J = 1, M
         S = A1 + B1*W2(J,2)
         W2(J,1) = -C1*W2(J,3)*K1(X,S)
  440 CONTINUE
      W2(1,1) = 1.0D0 + W2(1,1)
C     FIRST EQUATIONS
      IF ( .NOT. KCS) GO TO 520
      IF (GEVOD .AND. .NOT. GEV) GO TO 480
      DO 460 J = 1, M1
         MMJP1 = M - J + 1
         W1(1,J) = W2(J,1) + W2(MMJP1,1)
  460 CONTINUE
  480 IF (GEVOD .AND. GEV) GO TO 580
      DO 500 J = 1, M2
         MMJP1 = M - J + 1
         JDKS = ICM2 + J - 1
         W1(1,JDKS) = W2(J,1) - W2(MMJP1,1)
  500 CONTINUE
      GO TO 580
  520 DO 540 J = 1, M
         W1(1,J) = W2(J,1)
  540 CONTINUE
      X = A
      DO 560 J = 1, M
         S = A1 + B1*W2(J,2)
         W1(M,J) = -C1*W2(J,3)*K2(X,S)
  560 CONTINUE
      W1(M,M) = 1.0D0 + W1(M,M)
C     REMAINING EQUATIONS
  580 NR = 2.0D0*NR
      IF (N2P1.LT.2) GO TO 900
      DO 880 I = 2, N2P1
C        EVALUATE INTEGRALS
         T0 = 1.0D0
         T1 = W2(I,2)
         X2 = 2.0D0*T1
         WD(1) = T1 + 1.0D0
         V0 = 0.0D0
         V1 = 0.5D0*WD(1)
         SIG = -1.0D0
         I2 = NMAX + 2
         JODD = .TRUE.
         DO 640 J = 1, N1
            T2 = X2*T1 - T0
            V2 = 0.5D0*(T2+SIG)/DBLE(J+1)
            IF ( .NOT. JODD) GO TO 600
            WD(I2) = V2 - V0
            I2 = I2 + 1
            GO TO 620
  600       IDKS = I2 - NMAX - 1
            WD(IDKS) = V2 - V0
  620       SIG = -SIG
            T0 = T1
            T1 = T2
            V0 = V1
            V1 = V2
            JODD = .NOT. JODD
  640    CONTINUE
         I2 = I2 - 1
         IF ( .NOT. JODD) GO TO 660
         WD(N2P1) = 0.5D0*WD(N2P1)
         GO TO 680
  660    WD(I2) = 0.5D0*WD(I2)
C        FORM WEIGHTS
  680    DO 700 J = 1, N2P1
            V0 = NR*C06DBF(W2(J,2),WD,N2P1,2)
            V1 = NR*C06DBF(W2(J,2),WD(NMAX+2),I2-NMAX-1,3)
            W2(J,4) = V0 + V1
            N1MJP2 = N1 - J + 2
            W2(N1MJP2,4) = V0 - V1
  700    CONTINUE
         W2(1,4) = 0.5D0*W2(1,4)
         W2(M,4) = 0.5D0*W2(M,4)
C        FORM ITH ROW
         X = A1 + B1*W2(I,2)
         DO 720 J = 1, M
            S = A1 + B1*W2(J,2)
            W2(J,1) = -C1*(W2(J,4)*K1(X,S)+(W2(J,3)-W2(J,4))*K2(X,S))
  720    CONTINUE
         W2(I,1) = 1.0D0 + W2(I,1)
C        FORM ITH EQUATIONS
         I1 = I
         IF ( .NOT. KCS) GO TO 800
         IF (GEVOD .AND. .NOT. GEV) GO TO 760
         DO 740 J = 1, M1
            MMJP1 = M - J + 1
            W1(I1,J) = W2(J,1) + W2(MMJP1,1)
  740    CONTINUE
  760    IF ((GEVOD .AND. GEV) .OR. (I1.GT.M2)) GO TO 860
         DO 780 J = 1, M2
            MMJP1 = M - J + 1
            JDKS = ICM2 + J - 1
            W1(I1,JDKS) = W2(J,1) - W2(MMJP1,1)
  780    CONTINUE
         GO TO 860
  800    DO 820 J = 1, M
            W1(I1,J) = W2(J,1)
  820    CONTINUE
         I1 = N1 - I + 2
         X = A1 + B1*W2(I1,2)
         DO 840 J = 1, M
            S = A1 + B1*W2(J,2)
            J2 = N1 - J + 2
            W1(I1,J) = -C1*((W2(J2,3)-W2(J2,4))*K1(X,S)+W2(J2,4)*K2(X,S)
     *                 )
  840    CONTINUE
         W1(I1,I1) = 1.0D0 + W1(I1,I1)
  860    CONTINUE
  880 CONTINUE
C     END I LOOP
C     SOLVE EQUATIONS
      IF ((KCS .AND. GEVOD) .AND. .NOT. GEV) GO TO 920
  900 I = 1
      CALL F04AAF(W1,NMAX,W1(1,IGI1),NMAX,M1,1,W1(1,IGI1),NMAX,F,I)
      IF (I.NE.0) GO TO 1240
  920 IF (( .NOT. KCS) .OR. (GEVOD .AND. GEV)) GO TO 940
      I = 1
      CALL F04AAF(W1(1,ICM2),NMAX,W1(1,IGI2),NMAX,M2,1,W1(1,IGI2)
     *            ,NMAX,F,I)
      IF (I.NE.0) GO TO 1240
  940 IF (KCS) GO TO 980
      DO 960 I = 1, M
         F(I) = W1(I,IGI1)
  960 CONTINUE
      GO TO 1140
  980 IF (GEVOD) GO TO 1060
      DO 1040 I = 1, N2P1
         IF (I.LE.M2) GO TO 1000
         F(I) = W1(I,IGI1)
         GO TO 1020
 1000    F(I) = 0.5D0*(W1(I,IGI1)+W1(I,IGI2))
         N1MIP2 = N1 - I + 2
         F(N1MIP2) = 0.5D0*(W1(I,IGI1)-W1(I,IGI2))
 1020    CONTINUE
 1040 CONTINUE
      GO TO 1140
 1060 IF (GEV) GO TO 1100
      DO 1080 I1 = 1, M2
         F(I1) = 0.5D0*W1(I1,IGI2)
 1080 CONTINUE
      GO TO 1140
 1100 DO 1120 I = 1, N2P1
         F(I) = 0.5D0*W1(I,IGI1)
 1120 CONTINUE
C     FIND CHEBYSHEV COEFFICIENTS OF SOLUTION
 1140 IF ( .NOT. (KCS .AND. GEVOD)) F(N) = F(N)*0.5D0
      DO 1220 I = 1, N
         IF (KCS .AND. GEVOD) GO TO 1160
         X2 = COS(DBLE(I-1)*P)
         GO TO 1200
 1160    IF (GEV) GO TO 1180
         X2 = COS((2.0D0*DBLE(I-1)+1.0D0)*P)
         GO TO 1200
 1180    X2 = COS(2.0D0*DBLE(I-1)*P)
 1200    C(I) = NR*C06DBF(X2,F,N,1)
         IF (KCS .AND. GEVOD) C(I) = 2.0D0*C(I)
 1220 CONTINUE
      IF (KCS .AND. GEVOD) GO TO 1260
      C(N) = 0.5D0*C(N)
      F(N) = 2.0D0*F(N)
C     END GENERAL CASE
      GO TO 1260
 1240 IERROR = 2
 1260 IF (IERROR.EQ.0) IFAIL = 0
      IF (IERROR.NE.0) IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
