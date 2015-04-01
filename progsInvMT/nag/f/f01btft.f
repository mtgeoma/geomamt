      SUBROUTINE F01BTF(N,A,IA,P,DP,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1672 (JUN 1995).
C
C     DECOMPOSES A REAL MATRIX INTO A PRODUCT OF TRIANGULAR FACTORS
C     BY GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING.
C     THE ROUTINE OPERATES ON BLOCKS OF CONSECUTIVE COLUMNS
C     FOR EFFICIENCY ON A PAGED MACHINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BTF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DP
      INTEGER           IA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), P(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BASE, HALF, ONE, THRESH, X, Y, ZERO
      INTEGER           I, J, J1, K, K1, K2, L, LACT, LBLK
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF, X02BHF, X02CAZ
      EXTERNAL          X02AJF, P01ABF, X02BHF, X02CAZ
C     .. External Subroutines ..
      EXTERNAL          DGER, DSCAL, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              ZERO/0.0D0/, HALF/0.5D0/, ONE/1.0D0/
C     .. Executable Statements ..
C
      IF (N.LE.0 .OR. IA.LT.N) GO TO 260
C
C     THRESH IS THE THRESHOLD ON PIVOT SIZE RELATIVE TO THE ORIGINAL
C     INFINITY NORM OF ITS ROW.
C
      BASE = X02BHF()
      THRESH = BASE*X02AJF()*HALF
C
C     LBLK IS THE NUMBER OF COLUMNS THAT ARE TREATED TOGETHER.
C     LBLK IS CHOSEN SO THAT THERE IS ROOM IN THE DESIRED ACTIVE SET
C     FOR LBLK CONTIGUOUS COLUMNS OF A, ANOTHER COLUMN ON ITS OWN,
C     AND THE ARRAY P.
C
      LACT = X02CAZ(X)
      LBLK = N
      IF (LACT.NE.0) LBLK = MAX(1,LACT/IA-2)
C
C     FIND RECIPROCALS OF THE INFINITY NORMS OF THE ROWS OF A.
C
      DO 20 I = 1, N
         P(I) = ZERO
   20 CONTINUE
      DO 60 J = 1, N
         DO 40 I = 1, N
            P(I) = MAX(ABS(A(I,J)),P(I))
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, N
         IF (P(I).LE.ZERO) GO TO 220
         P(I) = ONE/P(I)
   80 CONTINUE
      DP = 1.0D0
C
C     MAIN LOOP IN WHICH A BLOCK OF COLUMNS IS TREATED.
C     THE ACTIVE BLOCK CONSISTS OF COLUMNS K1 TO K2.
C     THE PREVIOUS COLUMNS SHOULD BE PAGED IN AT MOST ONCE DURING
C     EXECUTION OF THIS LOOP.
C
      DO 200 K1 = 1, N, LBLK
         K2 = MIN(K1+LBLK-1,N)
         DO 180 K = 1, K2
            IF (K.LT.K1) GO TO 140
C           K-TH COLUMN IS IN ACTIVE BLOCK - CHOOSE PIVOT
            X = ZERO
            DO 100 I = K, N
               Y = ABS(A(I,K))*P(I)
               IF (Y.LE.X) GO TO 100
               X = Y
               L = I
  100       CONTINUE
            IF (X.LT.THRESH) GO TO 240
            IF (L.EQ.K) GO TO 120
C           PERFORM INTERCHANGE
            Y = A(K,K)
            A(K,K) = A(L,K)
            A(L,K) = Y
            P(L) = P(K)
            DP = -DP
  120       P(K) = L
            J1 = K + 1
            GO TO 160
C           K-TH COLUMN PRECEDES ACTIVE BLOCK - OBTAIN K-TH PIVOT INDEX
  140       L = P(K) + HALF
            J1 = K1
C           APPLY ELIMINATION OPERATIONS TO COLUMNS J1 TO K2 OF THE
C           ACTIVE BLOCK
  160       IF (K.LT.N) THEN
               IF (L.NE.K) CALL DSWAP(K2-J1+1,A(K,J1),IA,A(L,J1),IA)
               CALL DSCAL(K2-J1+1,ONE/A(K,K),A(K,J1),IA)
               CALL DGER(N-K,K2-J1+1,-ONE,A(K+1,K),1,A(K,J1),IA,
     *                   A(K+1,J1),IA)
            END IF
  180    CONTINUE
  200 CONTINUE
C
C     SUCCESSFUL EXIT
C
      IFAIL = 0
      GO TO 280
C
C     ERROR EXITS
C
  220 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      GO TO 280
  240 IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      GO TO 280
  260 IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
  280 RETURN
      END
