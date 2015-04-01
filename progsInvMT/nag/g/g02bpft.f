      SUBROUTINE G02BPF(N,M,X,IX,MISS,XMISS,ITYPE,RR,IRR,NCASES,INCASE,
     *                  KWORKA,KWORKB,KWORKC,WORK1,WORK2,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BPF
C     WRITTEN 30. 7.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES KENDALL AND/OR SPEARMAN RANK CORRELATION
C     COEFFICIENTS
C     FOR A SET OF DATA IN THE ARRAY X, OMITTING COMPLETELY ANY
C     CASES WITH A MISSING OBSERVATION FOR ANY VARIABLE -- THE
C     ARRAY X IS OVERWRITTEN BY THE ROUTINE
C
C
C     USES NAG ERROR ROUTINE P01AAF
C     NAG LIBRARY ROUTINE    X02BEF
C     AND AUXILIARY ROUTINES G02BNZ
C     G02BNX
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE-DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BPF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IRR, ITYPE, IX, M, N, NCASES
C     .. Array Arguments ..
      DOUBLE PRECISION  RR(IRR,M), WORK1(M), WORK2(M), X(IX,M), XMISS(M)
      INTEGER           INCASE(N), KWORKA(N), KWORKB(N), KWORKC(N),
     *                  MISS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, C, D, U
      INTEGER           I, IC, IERROR, IP, IV, IVP, J, JC, JV, K, NM
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF, X02BEF
      EXTERNAL          P01ABF, X02BEF
C     .. External Subroutines ..
      EXTERNAL          G02BNX, G02BNZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
      ACC = 0.1D0**(X02BEF(ACC)-2)
      IERROR = 0
      IF (ITYPE.LT.-1 .OR. ITYPE.GT.1) IERROR = 4
      IF (IX.LT.N .OR. IRR.LT.M) IERROR = 3
      IF (M.LT.2) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 760
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
      IF (NM.LE.0) GO TO 160
      IV = 0
      DO 140 I = 1, N
         INCASE(I) = 0
         DO 80 J = 1, NM
            JV = MISS(J)
            IF (ABS(X(I,JV)-XMISS(J)).LE.ABS(ACC*XMISS(J))) GO TO 100
   80    CONTINUE
         INCASE(I) = I
         IV = IV + 1
         KWORKC(IV) = I
         GO TO 140
  100    DO 120 J = 1, M
            X(I,J) = 0.0D0
  120    CONTINUE
  140 CONTINUE
      NCASES = IV
      IF (NCASES.LT.2) GO TO 740
      GO TO 200
  160 DO 180 I = 1, N
         INCASE(I) = I
         KWORKC(I) = I
  180 CONTINUE
      NCASES = N
  200 DO 260 IV = 1, M
         DO 220 I = 1, N
            KWORKA(I) = INCASE(I)
  220    CONTINUE
         DO 240 I = 1, NCASES
            KWORKB(I) = KWORKC(I)
  240    CONTINUE
         CALL G02BNZ(X(1,IV),N,NCASES,KWORKA,KWORKB)
         CALL G02BNX(X(1,IV),N,NCASES,KWORKA,KWORKB(1),WORK1(IV)
     *               ,WORK2(IV))
  260 CONTINUE
      NM = M - 1
      K = NCASES - 1
      IF (ITYPE) 280, 280, 420
  280 DO 400 IV = 1, NM
         IVP = IV + 1
         DO 380 JV = IVP, M
            U = 0.0D0
            DO 320 I = 1, K
               IC = KWORKC(I)
               IP = I + 1
               DO 300 J = IP, NCASES
                  JC = KWORKC(J)
                  C = X(IC,IV) - X(JC,IV)
                  D = X(IC,JV) - X(JC,JV)
                  IF (C.NE.0.0D0 .AND. D.NE.0.0D0) U = U + SIGN(1.0D0,C)
     *                *SIGN(1.0D0,D)
  300          CONTINUE
  320       CONTINUE
            D = WORK1(IV)*WORK1(JV)
            IF (D) 340, 340, 360
  340       RR(JV,IV) = 0.0D0
            GO TO 380
  360       RR(JV,IV) = U/SQRT(D)
  380    CONTINUE
  400 CONTINUE
  420 IF (ITYPE) 560, 440, 440
  440 DO 540 IV = 2, M
         IVP = IV - 1
         DO 520 JV = 1, IVP
            U = 0.0D0
            DO 460 I = 1, NCASES
               IC = KWORKC(I)
               D = X(IC,IV) - X(IC,JV)
               U = U + D*D
  460       CONTINUE
            D = WORK2(IV)*WORK2(JV)
            IF (D) 480, 480, 500
  480       RR(JV,IV) = 0.0D0
            GO TO 520
  500       RR(JV,IV) = 0.5D0*(WORK2(IV)+WORK2(JV)-6.0D0*U)/SQRT(D)
  520    CONTINUE
  540 CONTINUE
      IF (ITYPE) 560, 680, 620
  560 DO 600 IV = 2, M
         IVP = IV - 1
         DO 580 JV = 1, IVP
            RR(JV,IV) = RR(IV,JV)
  580    CONTINUE
  600 CONTINUE
      GO TO 680
  620 DO 660 IV = 1, NM
         IVP = IV + 1
         DO 640 JV = IVP, M
            RR(JV,IV) = RR(IV,JV)
  640    CONTINUE
  660 CONTINUE
  680 DO 700 I = 1, M
         RR(I,I) = 1.0D0
  700 CONTINUE
      DO 720 I = 1, N
         IF (INCASE(I).GT.0) INCASE(I) = 1
  720 CONTINUE
      IFAIL = 0
      RETURN
  740 IERROR = 5
  760 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
