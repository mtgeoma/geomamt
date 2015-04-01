      SUBROUTINE G02BSF(N,M,X,IX,MISS,XMISS,ITYPE,RR,IRR,NCASES,COUNT,
     *                  IC,KWORKA,KWORKB,KWORKC,KWORKD,WORK1,WORK2,
     *                  IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BSF
C     WRITTEN 16. 8.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES KENDALL AND/OR SPEARMAN RANK CORRELATION
C     COEFFICIENTS
C     FOR A SET OF DATA IN THE ARRAY X, OMITTING CASES WITH MISSING
C     VALUES FROM ONLY THOSE CALCULATIONS INVOLVING THE VARIABLES
C     FOR WHICH THE VALUES ARE MISSING--THE ARRAY X IS PRESERVED.
C
C     USES NAG ERROR ROUTINE P01AAF
C     NAG LIBRARY ROUTINE    X02BEF
C     AND AUXILIARY ROUTINES G02BNZ
C     G02BNY
C     G02BRZ
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE-DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BSF')
C     .. Scalar Arguments ..
      INTEGER           IC, IFAIL, IRR, ITYPE, IX, M, N, NCASES
C     .. Array Arguments ..
      DOUBLE PRECISION  COUNT(IC,M), RR(IRR,M), WORK1(N), WORK2(N),
     *                  X(IX,M), XMISS(M)
      INTEGER           KWORKA(N), KWORKB(N), KWORKC(N), KWORKD(N),
     *                  MISS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, FK
      INTEGER           I, IERROR, IP, IV, JV, K, KD, KV, M1
C     .. Local Arrays ..
      DOUBLE PRECISION  CORR(4)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF, X02BEF
      EXTERNAL          P01ABF, X02BEF
C     .. External Subroutines ..
      EXTERNAL          G02BNY, G02BNZ, G02BRZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      ACC = 0.1D0**(X02BEF(ACC)-2)
      IERROR = 0
      IF (ITYPE.LT.-1 .OR. ITYPE.GT.1) IERROR = 4
      IF (IX.LT.N .OR. IRR.LT.M .OR. IC.LT.M) IERROR = 3
      IF (M.LT.2) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 500
   20 NCASES = N
      IF (MISS(1)) 80, 80, 40
   40 KD = 0
      DO 60 I = 1, N
         KWORKC(I) = 0
         IF (ABS(X(I,1)-XMISS(1)).LE.ABS(ACC*XMISS(1))) GO TO 60
         KD = KD + 1
         KWORKC(I) = I
   60 CONTINUE
      GO TO 120
   80 KD = N
      DO 100 I = 1, N
         KWORKC(I) = I
  100 CONTINUE
  120 COUNT(1,1) = DBLE(KD)
      RR(1,1) = 1.0D0
      M1 = M - 1
      DO 460 IV = 1, M1
         DO 440 KV = IV, M1
            JV = M + IV - KV
            K = 0
            IF (KV.EQ.M1) GO TO 220
            IF (MISS(JV)) 180, 180, 140
  140       DO 160 I = 1, N
               KWORKA(I) = 0
               IF (ABS(X(I,JV)-XMISS(JV)).LE.ABS(ACC*XMISS(JV))
     *             .OR. KWORKC(I).EQ.0) GO TO 160
               K = K + 1
               KWORKA(I) = I
               KWORKB(K) = I
               KWORKD(K) = I
  160       CONTINUE
            GO TO 340
  180       DO 200 I = 1, N
               KWORKA(I) = KWORKC(I)
               IF (KWORKA(I).EQ.0) GO TO 200
               K = K + 1
               KWORKB(K) = I
               KWORKD(K) = I
  200       CONTINUE
            GO TO 340
  220       IF (MISS(JV)) 300, 300, 240
  240       KD = 0
            DO 280 I = 1, N
               IP = KWORKC(I)
               KWORKC(I) = 0
               KWORKA(I) = 0
               IF (ABS(X(I,JV)-XMISS(JV)).LE.ABS(ACC*XMISS(JV)))
     *             GO TO 260
               KWORKC(I) = I
               KD = KD + 1
  260          IF (KWORKC(I).LE.0 .OR. IP.LE.0) GO TO 280
               K = K + 1
               KWORKA(I) = I
               KWORKB(K) = I
               KWORKD(K) = I
  280       CONTINUE
            GO TO 340
  300       DO 320 I = 1, N
               KWORKA(I) = KWORKC(I)
               KWORKC(I) = I
               IF (KWORKA(I).EQ.0) GO TO 320
               K = K + 1
               KWORKB(K) = I
               KWORKD(K) = I
  320       CONTINUE
            KD = N
  340       IF (K.LT.2) GO TO 400
            CALL G02BNZ(X(1,IV),N,K,KWORKA,KWORKB)
            CALL G02BNY(X(1,IV),N,K,KWORKA,KWORKB(1),WORK1,CORR(1)
     *                  ,CORR(2))
            DO 360 I = 1, N
               KWORKA(I) = 0
  360       CONTINUE
            DO 380 I = 1, K
               IP = KWORKD(I)
               KWORKB(I) = KWORKD(I)
               KWORKA(IP) = IP
  380       CONTINUE
            CALL G02BNZ(X(1,JV),N,K,KWORKA,KWORKB)
            CALL G02BNY(X(1,JV),N,K,KWORKA,KWORKB(1),WORK2,CORR(3)
     *                  ,CORR(4))
            CALL G02BRZ(N,K,WORK1,WORK2,KWORKD,CORR,ITYPE,RR(JV,IV)
     *                  ,RR(IV,JV))
            GO TO 420
  400       RR(IV,JV) = 0.0D0
            RR(JV,IV) = 0.0D0
  420       FK = DBLE(K)
            COUNT(IV,JV) = FK
            COUNT(JV,IV) = FK
            IF (K.LT.NCASES) NCASES = K
  440    CONTINUE
         K = IV + 1
         RR(K,K) = 1.0D0
         COUNT(K,K) = DBLE(KD)
  460 CONTINUE
      IF (NCASES.LT.2) GO TO 480
      IFAIL = 0
      RETURN
  480 IERROR = 5
  500 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
