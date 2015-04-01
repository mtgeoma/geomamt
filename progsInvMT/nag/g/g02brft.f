      SUBROUTINE G02BRF(N,M,X,IX,MISS,XMISS,ITYPE,RR,IRR,NCASES,INCASE,
     *                  KWORKA,KWORKB,KWORKC,WORK1,WORK2,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BRF
C     WRITTEN 15. 8.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES KENDALL AND/OR SPEARMAN RANK CORRELATION
C     COEFFICIENTS
C     FOR A SET OF DATA IN THE ARRAY X, OMITTING COMPLETELY ANY
C     CASES WITH A MISSING OBSERVATION FOR ANY VARIABLE -- THE
C     ARRAY X IS PRESERVED.
C
C     USES NAG ERROR ROUTINE P01AAF
C     NAG LIBRARY ROUTINE    X02BEF
C     AND AUXILIARY ROUTINES G02BNZ
C     G02BNY
C     G02BRZ
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BRF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IRR, ITYPE, IX, M, N, NCASES
C     .. Array Arguments ..
      DOUBLE PRECISION  RR(IRR,M), WORK1(N), WORK2(N), X(IX,M), XMISS(M)
      INTEGER           INCASE(N), KWORKA(N), KWORKB(N), KWORKC(N),
     *                  MISS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC
      INTEGER           I, IERROR, IV, J, JV, K, KV, NM
C     .. Local Arrays ..
      DOUBLE PRECISION  CORR(4)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF, X02BEF
      EXTERNAL          P01ABF, X02BEF
C     .. External Subroutines ..
      EXTERNAL          G02BNY, G02BNZ, G02BRZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      ACC = 0.1D0**(X02BEF(ACC)-2)
      IERROR = 0
      IF (ITYPE.LT.-1 .OR. ITYPE.GT.1) IERROR = 4
      IF (IX.LT.N .OR. IRR.LT.M) IERROR = 3
      IF (M.LT.2) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 320
C
C     RE-ARRANGE MISSING VALUE INFORMATION FOR EFFICIENCY
C
   20 NM = 0
      DO 60 I = 1, M
         IF (MISS(I)) 60, 60, 40
   40    NM = NM + 1
         MISS(NM) = I
         XMISS(NM) = XMISS(I)
   60 CONTINUE
      IF (NM.LE.0) GO TO 120
      K = 0
      DO 100 I = 1, N
         KWORKA(I) = 0
         INCASE(I) = 0
         DO 80 J = 1, NM
            JV = MISS(J)
            IF (ABS(X(I,JV)-XMISS(J)).LE.ABS(ACC*XMISS(J))) GO TO 100
   80    CONTINUE
         KWORKA(I) = I
         INCASE(I) = I
         K = K + 1
         KWORKB(K) = I
         KWORKC(K) = I
  100 CONTINUE
      NCASES = K
      IF (NCASES.LT.2) GO TO 300
      GO TO 160
  120 DO 140 I = 1, N
         KWORKA(I) = I
         KWORKB(I) = I
         INCASE(I) = I
         KWORKC(I) = I
  140 CONTINUE
      NCASES = N
  160 CALL G02BNZ(X(1,1),N,NCASES,KWORKA,KWORKB)
      CALL G02BNY(X(1,1),N,NCASES,KWORKA,KWORKB(1),WORK1,CORR(1),CORR(2)
     *            )
      NM = M - 1
      DO 260 IV = 1, NM
         DO 220 KV = IV, NM
            JV = M + IV - KV
            DO 180 I = 1, N
               KWORKA(I) = INCASE(I)
  180       CONTINUE
            DO 200 I = 1, NCASES
               KWORKB(I) = KWORKC(I)
  200       CONTINUE
            CALL G02BNZ(X(1,JV),N,NCASES,KWORKA,KWORKB)
            CALL G02BNY(X(1,JV),N,NCASES,KWORKA,KWORKB(1),WORK2,CORR(3)
     *                  ,CORR(4))
            CALL G02BRZ(N,NCASES,WORK1,WORK2,KWORKC,CORR,ITYPE,RR(JV,IV)
     *                  ,RR(IV,JV))
  220    CONTINUE
         RR(IV,IV) = 1.0D0
         IF (IV.EQ.NM) GO TO 260
         DO 240 I = 1, N
            IF (INCASE(I).GT.0) WORK1(I) = WORK2(I)
  240    CONTINUE
         CORR(1) = CORR(3)
         CORR(2) = CORR(4)
  260 CONTINUE
      RR(M,M) = 1.0D0
      DO 280 I = 1, N
         IF (INCASE(I).GT.0) INCASE(I) = 1
  280 CONTINUE
      IFAIL = 0
      RETURN
  300 IERROR = 5
  320 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
