      SUBROUTINE G02BQF(N,M,X,IX,ITYPE,RR,IRR,KWORKA,KWORKB,WORK1,WORK2,
     *                  IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BQF
C     WRITTEN 30. 7.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES KENDALL AND/OR SPEARMAN RANK CORRELATION
C     COEFFICIENTS
C     FOR A SET OF DATA IN THE ARRAY X -- THE ARRAY X IS PRESERVED.
C
C     USES NAG ERROR ROUTINE P01AAF
C     AND AUXILIARY ROUTINES G02BNZ
C     G02BNY
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE-DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BQF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IRR, ITYPE, IX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RR(IRR,M), WORK1(N), WORK2(N), X(IX,M)
      INTEGER           KWORKA(N), KWORKB(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DK, DS, U
      INTEGER           I, IERROR, IP, IV, J, JV, KV, M1, N1
C     .. Local Arrays ..
      DOUBLE PRECISION  CORR(4)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02BNY, G02BNZ
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN, SQRT
C     .. Executable Statements ..
      IERROR = 0
      IF (ITYPE.LT.-1 .OR. ITYPE.GT.1) IERROR = 4
      IF (IX.LT.N .OR. IRR.LT.M) IERROR = 3
      IF (M.LT.2) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 460
   20 M1 = M - 1
      N1 = N - 1
      DO 40 I = 1, N
         KWORKA(I) = I
         KWORKB(I) = I
   40 CONTINUE
      CALL G02BNZ(X(1,1),N,N,KWORKA,KWORKB)
      CALL G02BNY(X(1,1),N,N,KWORKA,KWORKB(1),WORK1,CORR(1),CORR(2))
      DO 440 IV = 1, M1
         DO 400 KV = IV, M1
            JV = M + IV - KV
            DO 60 I = 1, N
               KWORKA(I) = I
               KWORKB(I) = I
   60       CONTINUE
            CALL G02BNZ(X(1,JV),N,N,KWORKA,KWORKB)
            CALL G02BNY(X(1,JV),N,N,KWORKA,KWORKB(1),WORK2,CORR(3)
     *                  ,CORR(4))
            U = 0.0D0
            D = 0.0D0
            DO 100 I = 1, N1
               IP = I + 1
               DO 80 J = IP, N
                  DK = WORK1(I) - WORK1(J)
                  DS = WORK2(I) - WORK2(J)
                  IF (DK.NE.0.0D0 .AND. DS.NE.0.0D0) U = U +
     *                SIGN(1.0D0,DK)*SIGN(1.0D0,DS)
   80          CONTINUE
               DK = WORK1(I) - WORK2(I)
               D = D + DK*DK
  100       CONTINUE
            DK = WORK1(N) - WORK2(N)
            D = D + DK*DK
            IF (ITYPE) 120, 200, 320
  120       DK = CORR(1)*CORR(3)
            IF (DK) 140, 140, 160
  140       RR(JV,IV) = 0.0D0
            GO TO 180
  160       RR(JV,IV) = U/SQRT(DK)
  180       RR(IV,JV) = RR(JV,IV)
            GO TO 400
  200       DK = CORR(1)*CORR(3)
            IF (DK) 220, 220, 240
  220       RR(JV,IV) = 0.0D0
            GO TO 260
  240       RR(JV,IV) = U/SQRT(DK)
  260       DS = CORR(2)*CORR(4)
            IF (DS) 280, 280, 300
  280       RR(IV,JV) = 0.0D0
            GO TO 400
  300       RR(IV,JV) = 0.5D0*(CORR(2)+CORR(4)-6.0D0*D)/SQRT(DS)
            GO TO 400
  320       DS = CORR(2)*CORR(4)
            IF (DS) 340, 340, 360
  340       RR(JV,IV) = 0.0D0
            GO TO 380
  360       RR(JV,IV) = 0.5D0*(CORR(2)+CORR(4)-6.0D0*D)/SQRT(DS)
  380       RR(IV,JV) = RR(JV,IV)
  400    CONTINUE
         RR(IV,IV) = 1.0D0
         IF (IV.EQ.M1) GO TO 440
         DO 420 I = 1, N
            WORK1(I) = WORK2(I)
  420    CONTINUE
         CORR(1) = CORR(3)
         CORR(2) = CORR(4)
  440 CONTINUE
      RR(M,M) = 1.0D0
      IFAIL = 0
      RETURN
  460 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
