      SUBROUTINE F01AVF(N,IB,AR,IAR,AI,IAI,LOW,LHI,D)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMBALANCE
C     REDUCE THE NORM OF A COMPLEX MATRIX AR,AI(N,N) BY
C     EXACT DIAGONAL SIMILARITY TRANSFORMATIONS STORED IN D(N).
C     DECEMBER 1ST.,1971
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IB, LHI, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), D(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B2, C, F, G, R, S
      INTEGER           I, J, JJ, K, L
      LOGICAL           NOCONV
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF
      EXTERNAL          A02ABF
C     .. External Subroutines ..
      EXTERNAL          F01AVZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      B2 = IB*IB
      L = 1
      K = N
   20 IF (K.LT.1) GO TO 100
C     SEARCH FOR ROWS ISOLATING AN EIGENVALUE AND PUSH THEM DOWN.
      J = K + 1
      DO 60 JJ = 1, K
         J = J - 1
         R = 0.0D0
         DO 40 I = 1, K
            IF (I.EQ.J) GO TO 40
            R = R + A02ABF(AR(J,I),AI(J,I))
   40    CONTINUE
         IF (R.EQ.0.0D0) GO TO 80
   60 CONTINUE
      GO TO 100
   80 CALL F01AVZ(K,AR,IAR,AI,IAI,D,K,L,N,J)
      K = K - 1
      GO TO 20
  100 IF (L.GT.K) GO TO 180
C     SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE AND PUSH THEM
C     LEFT.
      DO 140 J = L, K
         C = 0.0D0
         DO 120 I = L, K
            IF (I.EQ.J) GO TO 120
            C = C + A02ABF(AR(I,J),AI(I,J))
  120    CONTINUE
         IF (C.EQ.0.0D0) GO TO 160
  140 CONTINUE
      GO TO 180
  160 CALL F01AVZ(L,AR,IAR,AI,IAI,D,K,L,N,J)
      L = L + 1
      GO TO 100
C     NOW BALANCE THE SUB-MATRIX IN ROWS L THROUGH K.
  180 LOW = L
      LHI = K
      IF (L.GT.K) GO TO 220
      DO 200 I = L, K
         D(I) = 1.0D0
  200 CONTINUE
  220 NOCONV = .FALSE.
      IF (L.GT.K) GO TO 420
      DO 400 I = L, K
         C = 0.0D0
         R = 0.0D0
         DO 240 J = L, K
            IF (J.EQ.I) GO TO 240
            C = C + A02ABF(AR(J,I),AI(J,I))
            R = R + A02ABF(AR(I,J),AI(I,J))
  240    CONTINUE
         G = R/DBLE(IB)
         F = 1.0D0
         S = C + R
  260    IF (C.GE.G) GO TO 280
         F = F*DBLE(IB)
         C = C*B2
         GO TO 260
  280    G = R*DBLE(IB)
  300    IF (C.LT.G) GO TO 320
         F = F/DBLE(IB)
         C = C/B2
         GO TO 300
  320    IF (((C+R)/F).GE.(0.95D0*S)) GO TO 400
         G = 1.0D0/F
         D(I) = D(I)*F
         NOCONV = .TRUE.
         IF (L.GT.N) GO TO 360
         DO 340 J = L, N
            AR(I,J) = AR(I,J)*G
            AI(I,J) = AI(I,J)*G
  340    CONTINUE
  360    DO 380 J = 1, K
            AR(J,I) = AR(J,I)*F
            AI(J,I) = AI(J,I)*F
  380    CONTINUE
  400 CONTINUE
  420 IF (NOCONV) GO TO 220
      RETURN
      END
