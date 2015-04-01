      SUBROUTINE G08RBZ(Y,N,IRANK,NUNC,GAMMA,TOL,ZIN,ETA,VAPVEC,N1,ICO,
     *                  IWA,NIWA)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATES EXPECTED VALUES AND VARIANCE-COVARIANCE MATRIX
C     OF PARTICULAR FUNCTION OF ORDER STATISTICS FOR A GENERALISED
C     LOGISTIC DISTRIBUTION WHEN SOME OBSERVATIONS ARE RIGHT CENSORED.
C
C     PETTITT A.N.  APPROXIMATE METHODS USING RANKS FOR REGRESSION
C                   WITH CENSORED DATA.
C                   BIOMETRIKA, 70, PP 121-32.
C
C     ARGUMENTS :
C              Y - VECTOR OF OBSERVATIONS ( OR RANKS ).
C              N - SAMPLE SIZE.
C          IRANK - ARRAY USED TO STORE RANKS OF OBSERVATIONS.
C           NUNC - SPECIFIES NUMBER OF UNCENSORED OBSERVATIONS.
C          GAMMA - PARAMETER DEFINING THE GENERALIED LOGISTIC
C                  DISTRIBUTION
C                  GAMMA = ( 1 - LOGISTIC
C                          ( LIM - EXTREME VALUE
C                          ( APPROACHING 0
C            ZIN - ON EXIT CONTAINS EXPECTED VALUE OF FUNCTION OF ORDER
C                  STATISTICS IN ORDER CORRESPONDING TO RANKS OF
C                  OBSERVATIONS.
C            ETA - ON EXIT CONTAINS EXPECTED VALUE OF DERIVATIVE OF
C                  FUNCTION OF ORDER STATISTIC IN ORDER CORRESPONDING
C                  TO RANKS OF OBSERVATIONS.
C         VAPVEC - ON EXIT CONTAINS VARIANCE-COVARIANCE MATRIX OF
C                  FUNCTION OF ORDER STATISTICS IN ORDER CORRESPONDING
C                  TO RANKS OF OBSERVATIONS.
C             N1 - DIMENSION OF VECTOR NEEDED TO STORE UPPER TRIANGLE OF
C                  SQUARE SYMMETRIC MATRIX OF DIMENSION N.
C                  ( N1 .GE. N*(N+1)/2 ).
C            ICO - ON EXIT SPECIFIES THE NUMBER OF SETS
C                  OF TIES IN THE DATA
C            IWA - ON ENTRY, IWA CONTAINS AN INDEX TO THE ORIGINAL
C                  OBSERVATIONS IN Y AND SPECIFIES WHETHER OR NOT
C                  EACH Y(I) IS CENSORED. ON EXIT IT ALSO
C                  SPECIFIES THE NUMBER OF CENSORED OBSERVATIONS
C                  BETWEEN EACH PAIR OF UNCENSORED OBSERVATIONS AND
C                  THE LOCATION AND LENGTH OF SETS OF TIES.
C           NIWA - ON ENTRY, NIWA SPECIFIES THE LENGTH OF IWA.
C                  NIWA .GE. 4*N
C
C
C     ASSIGN CORRECT RANK TO THE CENSORED OBSERVATIONS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA, TOL
      INTEGER           ICO, N, N1, NIWA, NUNC
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), VAPVEC(N1), Y(N), ZIN(N)
      INTEGER           IRANK(N), IWA(NIWA)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, AN, AN1, AN2, GAM2, P1, P2, P3, SUM, SUM2,
     *                  ZI, ZJ
      INTEGER           I, ICI, ICJ, IFOUND, II, INC, INDI, INDJ, IRJ,
     *                  J, JC, JJ, K, M, N2, N3, NNSUM, NSUM
C     .. External Functions ..
      INTEGER           G01DCU
      EXTERNAL          G01DCU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      N2 = 2*N
      N3 = 3*N
      II = 0
      DO 60 I = 1, N
         IF (IWA(N+IWA(I)).EQ.0) THEN
            II = II + 1
         ELSE
            J = I
            IFOUND = 0
   20       IF (IFOUND.NE.0 .OR. J.EQ.N) GO TO 40
            J = J + 1
            IF (Y(J).GT.Y(I)) IFOUND = 1
            IF (Y(J).EQ.Y(I) .AND. IWA(N+IWA(J)).EQ.0) IFOUND = 2
            GO TO 20
   40       IF (IFOUND.EQ.2) THEN
               INDI = IWA(I)
               IWA(I) = IWA(J)
               IWA(J) = INDI
               II = II + 1
            END IF
         END IF
         IRANK(IWA(I)) = II
   60 CONTINUE
      I = 0
   80 I = I + 1
      IF (IWA(N+IWA(I)).EQ.1) GO TO 80
      JJ = I
      NNSUM = 0
      DO 140 I = 1, NUNC
         II = 0
  100    JJ = JJ + 1
         IF (JJ.EQ.N+1) GO TO 120
         IF (IWA(N+IWA(JJ)).EQ.0) GO TO 120
         II = II + 1
         GO TO 100
  120    IWA(N2+I) = II
         NNSUM = NNSUM + IWA(N2+I) + 1
  140 CONTINUE
C
C     AMEND RANKING IN CASE OF TIES AMONG UNCENSORED OBSERVATIONS
C
      II = 1
      JJ = 1
      ICO = 0
      J = 0
  160 IF (J.GE.NUNC-1) GO TO 260
      JC = 0
  180 J = J + 1
      JC = JC + 1
      IF (J.EQ.NUNC) GO TO 200
      JJ = JJ + IWA(N2+J) + 1
      IF (ABS(Y(II)-Y(JJ)).LE.TOL) GO TO 180
  200 IF (JC.EQ.1) THEN
         II = JJ
         GO TO 160
      ELSE
         JJ = II
         NSUM = 0
         DO 240 I = 1, JC - 1
            JJ = JJ + IWA(N2+J-JC+I) + 1
            IF (IWA(N2+J-JC+I).EQ.0) GO TO 240
            A = Y(JJ)
            IRJ = IRANK(IWA(JJ))
            INDJ = IWA(JJ)
            DO 220 K = 1, IWA(N2+J-JC+I)
               Y(JJ-K+1) = Y(JJ-K)
               IWA(JJ-K+1) = IWA(JJ-K)
               IRANK(IWA(JJ-K+1)) = IRJ
  220       CONTINUE
            Y(II+1) = A
            IWA(II+I) = INDJ
            IRANK(IWA(II+I)) = IRJ
            NSUM = NSUM + IWA(N2+J-JC+I)
            IWA(N2+J-JC+I) = 0
  240    CONTINUE
         IWA(N3+2*ICO+1) = II
         IWA(N3+2*ICO+2) = JC
         ICO = ICO + 1
         II = JJ + IWA(N2+J) + 1
         JJ = II
         IWA(N2+J) = IWA(N2+J) + NSUM
         GO TO 160
      END IF
C
C     CALCULATE SCORES FOR DIFFERENT VALUES OF GAMMA, FOR
C     GAMMA LESS THAN 0.0001 THE EXTREME VALUE SCORES ARE USED.
C
  260 IF (GAMMA.LE.0.0001D0) GO TO 340
      GAM2 = 2*GAMMA
      P1 = 1.0D0
      JJ = 0
      AN1 = NNSUM
      DO 320 J = 1, N
         IF (IWA(N+IWA(J)).EQ.0) THEN
            JJ = JJ + 1
            P1 = P1*AN1/(AN1+GAMMA)
            ZJ = P1
            ZIN(IWA(J)) = P1
            AN1 = AN1 - (IWA(N2+JJ)+1)
         ELSE
            ZIN(IWA(J)) = (1-P1)/GAMMA
         END IF
         P2 = 1.0D0
         P3 = P1
         II = 0
         INC = 0
         AN2 = NNSUM
         DO 300 I = 1, J
            IF (IRANK(IWA(I)).EQ.0 .OR. IRANK(IWA(J)).EQ.0) THEN
               M = G01DCU(IWA(I),IWA(J))
               VAPVEC(M) = 0.0D0
               GO TO 300
            END IF
            IF (IWA(N+IWA(J)).EQ.1 .AND. IRANK(IWA(I)).EQ.II)
     *          GO TO 280
            IF (IWA(N+IWA(J)).EQ.0 .AND. IWA(N+IWA(I)).EQ.1) GO TO 280
            IF (INC.EQ.1) THEN
               AN2 = AN2 - (IWA(N2+II)+1)
            END IF
            II = II + 1
            P2 = P2*AN2/(AN2+GAM2)
            P3 = P3*(AN2+GAMMA)/AN2
            ZI = ZIN(IWA(I))
  280       M = G01DCU(IWA(I),IWA(J))
            ICI = IWA(N+IWA(I))
            ICJ = IWA(N+IWA(J))
            VAPVEC(M) = (1/GAMMA+1-ICI)*(1/GAMMA+1-ICJ)*(P2*P3-ZI*ZJ)
            IF (IWA(N+IWA(I)).EQ.0) INC = 1
  300    CONTINUE
         IF (IWA(N+IWA(J)).EQ.0) THEN
            ETA(IWA(J)) = P1 - P2
         ELSE
            ETA(IWA(J)) = (P1-P2)/GAMMA
         END IF
  320 CONTINUE
      GO TO 400
C
C     EXTREME VALUE
C
  340 AN = NNSUM
      SUM = 0.0D0
      SUM2 = 0.0D0
      II = 0
      DO 380 I = 1, N
         IF (IWA(N+IWA(I)).EQ.0) THEN
            II = II + 1
            SUM = SUM + 1.0D0/AN
            SUM2 = SUM2 + 1.0D0/(AN*AN)
            AN = AN - (IWA(N2+II)+1)
         END IF
         ZIN(IWA(I)) = SUM
         ETA(IWA(I)) = SUM
         JJ = II - 1
         DO 360 J = I, N
            M = G01DCU(IWA(I),IWA(J))
            VAPVEC(M) = SUM2
  360    CONTINUE
  380 CONTINUE
  400 RETURN
      END
