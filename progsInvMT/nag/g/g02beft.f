      SUBROUTINE G02BEF(N,M,X,IX,MISS,XMISS,XBAR,STD,SSPZ,ISSPZ,RZ,IRZ,
     *                  NCASES,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BEF
C     WRITTEN 17. 7.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES MEANS AND STANDARD DEVIATIONS OF VARIABLES, SUMS OF
C     SQUARES AND CROSS-PRODUCTS ABOUT ZERO AND CORRELATION-LIKE
C     COEFFICIENTS FOR A SET OF DATA IN THE ARRAY X, OMITTING
C     COMPLETELY ANY CASES WITH A MISSING OBSERVATION FOR ANY
C     VARIABLE
C
C     USES NAG ERROR ROUTINE P01AAF
C     NAG LIBRARY ROUTINE    X02BEF
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE-DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BEF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IRZ, ISSPZ, IX, M, N, NCASES
C     .. Array Arguments ..
      DOUBLE PRECISION  RZ(IRZ,M), SSPZ(ISSPZ,M), STD(M), X(IX,M),
     *                  XBAR(M), XMISS(M)
      INTEGER           MISS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, S, XCASES, XM, XN
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
      IF (IX.LT.N .OR. ISSPZ.LT.M .OR. IRZ.LT.M) IERROR = 3
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
C     SOME MISSING VALUES ARE INVOLVED
C
   80 DO 160 I = 1, N
         DO 100 K = 1, NM
            J = MISS(K)
            IF (ABS(X(I,J)-XMISS(K)).LE.ABS(ACC*XMISS(K))) GO TO 160
  100    CONTINUE
         ISTART = I + 1
         DO 140 J = 1, M
            DO 120 K = J, M
               SSPZ(K,J) = X(I,J)*X(I,K)
  120       CONTINUE
            STD(J) = 0.0D0
            XBAR(J) = X(I,J)
  140    CONTINUE
         GO TO 300
  160 CONTINUE
      GO TO 580
C
C     NO MISSING VALUES INVOLVED IN FACT
C
  180 XCASES = DBLE(N)
      DO 280 J = 1, M
         S = 0.0D0
         DO 200 I = 1, N
            S = S + X(I,J)
  200    CONTINUE
         XM = S/XCASES
         XBAR(J) = XM
         S = 0.0D0
         DO 220 I = 1, N
            XN = X(I,J) - XM
            S = S + XN*XN
  220    CONTINUE
         STD(J) = S
         DO 260 K = J, M
            S = 0.0D0
            DO 240 I = 1, N
               S = S + X(I,J)*X(I,K)
  240       CONTINUE
            SSPZ(K,J) = S
  260    CONTINUE
  280 CONTINUE
      GO TO 400
C
C     SECOND AND SUBSEQUENT CASES WHEN MISSING VALUES ARE INVOLVED
C
  300 IF (ISTART.GT.N) GO TO 600
      XCASES = 1.0D0
      DO 380 I = ISTART, N
         DO 320 K = 1, NM
            J = MISS(K)
            IF (ABS(X(I,J)-XMISS(K)).LE.ABS(ACC*XMISS(K))) GO TO 380
  320    CONTINUE
         XCASES = XCASES + 1.0D0
         S = (XCASES-1.0D0)/XCASES
         DO 360 J = 1, M
            DO 340 K = J, M
               SSPZ(K,J) = SSPZ(K,J) + X(I,J)*X(I,K)
  340       CONTINUE
            XM = X(I,J) - XBAR(J)
            STD(J) = STD(J) + XM*XM*S
            XBAR(J) = XBAR(J) + XM/XCASES
  360    CONTINUE
  380 CONTINUE
      IF (XCASES.LT.1.5D0) GO TO 600
C
C     DO THIS WHETHER OR NOT MISSING VALUES ARE INVOLVED
C
  400 XN = XCASES - 1.0D0
      DO 560 J = 1, M
         IF (STD(J)) 420, 420, 440
  420    STD(J) = 0.0D0
         GO TO 460
  440    STD(J) = SQRT(STD(J)/XN)
  460    DO 540 K = 1, J
            S = SSPZ(J,J)*SSPZ(K,K)
            IF (S) 480, 480, 500
  480       RZ(J,K) = 0.0D0
            GO TO 520
  500       RZ(J,K) = SSPZ(J,K)/SQRT(S)
  520       RZ(K,J) = RZ(J,K)
            SSPZ(K,J) = SSPZ(J,K)
  540    CONTINUE
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
