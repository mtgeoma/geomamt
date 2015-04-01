      SUBROUTINE G02BBF(N,M,X,IX,MISS,XMISS,XBAR,STD,SSP,ISSP,R,IR,
     *                  NCASES,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BBF
C     WRITTEN 16. 7.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES MEANS, STANDARD DEVIATIONS, SUMS OF SQUARES AND
C     CROSS-PRODUCTS OF DEVIATIONS FROM MEANS, AND PEARSON PRODUCT-
C     MOMENT CORRELATION COEFFICIENTS FOR A SET OF DATA IN THE
C     ARRAY X, OMITTING COMPLETELY ANY CASES WITH A MISSING
C     OBSERVATION FOR ANY VARIABLE
C
C     USES NAG ERROR ROUTINE P01AAF
C     NAG LIBRARY ROUTINE    X02BEF
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH MAY BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IR, ISSP, IX, M, N, NCASES
C     .. Array Arguments ..
      DOUBLE PRECISION  R(IR,M), SSP(ISSP,M), STD(M), X(IX,M), XBAR(M),
     *                  XMISS(M)
      INTEGER           MISS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, S, XCASES
      INTEGER           I, IERROR, ISTART, J, K, NM
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF, X02BEF
      EXTERNAL          P01ABF, X02BEF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT, DBLE, INT
C     .. Executable Statements ..
      ACC = 0.1D0**(X02BEF(ACC)-2)
      IERROR = 0
      IF (IX.LT.N .OR. ISSP.LT.M .OR. IR.LT.M) IERROR = 3
      IF (M.LT.2) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 620
   20 XCASES = 0.0D0
C
C     RE-ARRANGE MISSING VALUE INFORMATION FOR EFFICIENCY
C
      NM = 0
      DO 60 I = 1, M
         IF (MISS(I)) 60, 60, 40
   40    NM = NM + 1
         MISS(NM) = I
         XMISS(NM) = XMISS(I)
   60 CONTINUE
      IF (NM) 180, 180, 80
C
C     SOME MISSING VALUES ARE SPECIFIED
C
   80 DO 160 I = 1, N
         DO 100 K = 1, NM
            J = MISS(K)
            IF (ABS(X(I,J)-XMISS(K)).LE.ABS(ACC*XMISS(K))) GO TO 160
  100    CONTINUE
         ISTART = I + 1
         DO 140 J = 1, M
            DO 120 K = J, M
               SSP(K,J) = 0.0D0
  120       CONTINUE
            XBAR(J) = X(I,J)
  140    CONTINUE
         GO TO 280
  160 CONTINUE
      GO TO 580
C
C     NO MISSING VALUES INVOLVED IN FACT
C
  180 XCASES = DBLE(N)
      DO 260 J = 1, M
         S = 0.0D0
         DO 200 I = 1, N
            S = S + X(I,J)
  200    CONTINUE
         XBAR(J) = S/XCASES
         DO 240 K = 1, J
            S = 0.0D0
            DO 220 I = 1, N
               S = S + (X(I,J)-XBAR(J))*(X(I,K)-XBAR(K))
  220       CONTINUE
            SSP(J,K) = S
  240    CONTINUE
  260 CONTINUE
      GO TO 380
C
C     SECOND AND SUBSEQUENT CASES WHEN MISSING VALUES ARE SPECIFIED
C
  280 IF (ISTART.GT.N) GO TO 600
      XCASES = 1.0D0
      DO 360 I = ISTART, N
         DO 300 K = 1, NM
            J = MISS(K)
            IF (ABS(X(I,J)-XMISS(K)).LE.ABS(ACC*XMISS(K))) GO TO 360
  300    CONTINUE
         XCASES = XCASES + 1.0D0
         S = (XCASES-1.0D0)/XCASES
         DO 340 J = 1, M
            DO 320 K = J, M
               SSP(K,J) = SSP(K,J) + (X(I,J)-XBAR(J))*(X(I,K)-XBAR(K))*S
  320       CONTINUE
            XBAR(J) = XBAR(J)*S + X(I,J)/XCASES
  340    CONTINUE
  360 CONTINUE
      IF (XCASES.LT.1.5D0) GO TO 600
C
C     DO THIS WHETHER OR NOT MISSING VALUES ARE INVOLVED
C
  380 DO 540 J = 1, M
         S = SSP(J,J)
         IF (S) 400, 400, 420
  400    STD(J) = 0.0D0
         GO TO 440
  420    STD(J) = SQRT(S)
  440    DO 520 K = 1, J
            S = STD(J)*STD(K)
            IF (S) 460, 460, 480
  460       R(J,K) = 0.0D0
            GO TO 500
  480       R(J,K) = SSP(J,K)/S
  500       R(K,J) = R(J,K)
            SSP(K,J) = SSP(J,K)
  520    CONTINUE
  540 CONTINUE
      S = SQRT(XCASES-1.0D0)
      DO 560 J = 1, M
         STD(J) = STD(J)/S
  560 CONTINUE
      NCASES = INT(XCASES+0.5D0)
      IFAIL = 0
      RETURN
  580 NCASES = 0
      IERROR = 4
      GO TO 620
  600 NCASES = 1
      IERROR = 5
  620 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
