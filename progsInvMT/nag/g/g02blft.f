      SUBROUTINE G02BLF(N,M,X,IX,MISS,XMISS,MISTYP,NVARS,KVAR,XBAR,STD,
     *                  SSPZ,ISSPZ,RZ,IRZ,NCASES,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BLF
C     WRITTEN 17.8.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES MEANS AND STANDARD DEVIATIONS OF VARIABLES, SUMS OF
C     SQUARES AND CROSS-PRODUCTS ABOUT ZERO AND CORRELATION-LIKE
C     COEFFICIENTS FOR A SET OF DATA IN SPECIFIED COLUMNS OF THE
C     ARRAY X, OMITTING COMPLETELY ANY CASES WITH A MISSING
C     OBSERVATION FOR ANY VARIABLE (OVER ALL M VARIABLES IF
C     MISTYP =1, OR OVER ONLY THE NVARS VARIABLES IN THE SELECTED
C     SUBSET IF MISTYP=0)
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
      PARAMETER         (SRNAME='G02BLF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IRZ, ISSPZ, IX, M, MISTYP, N, NCASES,
     *                  NVARS
C     .. Array Arguments ..
      DOUBLE PRECISION  RZ(IRZ,NVARS), SSPZ(ISSPZ,NVARS), STD(NVARS),
     *                  X(IX,M), XBAR(NVARS), XMISS(M)
      INTEGER           KVAR(NVARS), MISS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, S, XCASES, XM, XN
      INTEGER           I, IERROR, ISTART, J, JV, K, KV, NM
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
      IF (MISTYP.GT.1 .OR. MISTYP.LT.0) IERROR = 5
      IF (IX.LT.N .OR. ISSPZ.LT.NVARS .OR. IRZ.LT.NVARS) IERROR = 3
      IF (NVARS.LT.2 .OR. NVARS.GT.M) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 760
   20 DO 40 I = 1, NVARS
         IF (KVAR(I).LT.1 .OR. KVAR(I).GT.M) GO TO 700
   40 CONTINUE
      XCASES = 0.0D0
C
C     RE-ARRANGE MISSING VALUE INFORMATION FOR EFFICIENCY
C
      NM = 0
      IF (MISTYP) 60, 60, 120
   60 DO 100 JV = 1, NVARS
         I = KVAR(JV)
         IF (MISS(I)) 100, 100, 80
   80    NM = NM + 1
         MISS(NM) = I
         XMISS(NM) = XMISS(I)
  100 CONTINUE
      GO TO 180
  120 DO 160 I = 1, M
         IF (MISS(I)) 160, 160, 140
  140    NM = NM + 1
         MISS(NM) = I
         XMISS(NM) = XMISS(I)
  160 CONTINUE
  180 IF (NM) 300, 300, 200
C
C     SOME MISSING VALUES ARE INVOLVED
C
  200 DO 280 I = 1, N
         DO 220 K = 1, NM
            J = MISS(K)
            IF (ABS(X(I,J)-XMISS(K)).LE.ABS(ACC*XMISS(K))) GO TO 280
  220    CONTINUE
         ISTART = I + 1
         DO 260 J = 1, NVARS
            JV = KVAR(J)
            DO 240 K = J, NVARS
               KV = KVAR(K)
               SSPZ(K,J) = X(I,JV)*X(I,KV)
  240       CONTINUE
            STD(J) = 0.0D0
            XBAR(J) = X(I,JV)
  260    CONTINUE
         GO TO 420
  280 CONTINUE
      GO TO 720
C
C     NO MISSING VALUES INVOLVED IN FACT
C
  300 XCASES = DBLE(N)
      DO 400 J = 1, NVARS
         JV = KVAR(J)
         S = 0.0D0
         DO 320 I = 1, N
            S = S + X(I,JV)
  320    CONTINUE
         XM = S/XCASES
         XBAR(J) = XM
         S = 0.0D0
         DO 340 I = 1, N
            XN = X(I,JV) - XM
            S = S + XN*XN
  340    CONTINUE
         STD(J) = S
         DO 380 K = J, NVARS
            KV = KVAR(K)
            S = 0.0D0
            DO 360 I = 1, N
               S = S + X(I,JV)*X(I,KV)
  360       CONTINUE
            SSPZ(K,J) = S
  380    CONTINUE
  400 CONTINUE
      GO TO 520
C
C     SECOND AND SUBSEQUENT CASES WHEN MISSING VALUES ARE INVOLVED
C
  420 IF (ISTART.GT.N) GO TO 740
      XCASES = 1.0D0
      DO 500 I = ISTART, N
         DO 440 K = 1, NM
            J = MISS(K)
            IF (ABS(X(I,J)-XMISS(K)).LE.ABS(ACC*XMISS(K))) GO TO 500
  440    CONTINUE
         XCASES = XCASES + 1.0D0
         S = (XCASES-1.0D0)/XCASES
         DO 480 J = 1, NVARS
            JV = KVAR(J)
            DO 460 K = J, NVARS
               KV = KVAR(K)
               SSPZ(K,J) = SSPZ(K,J) + X(I,JV)*X(I,KV)
  460       CONTINUE
            XM = X(I,JV) - XBAR(J)
            STD(J) = STD(J) + XM*XM*S
            XBAR(J) = XBAR(J) + XM/XCASES
  480    CONTINUE
  500 CONTINUE
      IF (XCASES.LT.1.5D0) GO TO 740
C
C     DO THIS WHETHER OR NOT MISSING VALUES ARE INVOLVED
C
  520 XN = XCASES - 1.0D0
      DO 680 J = 1, NVARS
         IF (STD(J)) 540, 540, 560
  540    STD(J) = 0.0D0
         GO TO 580
  560    STD(J) = SQRT(STD(J)/XN)
  580    DO 660 K = 1, J
            S = SSPZ(J,J)*SSPZ(K,K)
            IF (S) 600, 600, 620
  600       RZ(J,K) = 0.0D0
            GO TO 640
  620       RZ(J,K) = SSPZ(J,K)/SQRT(S)
  640       RZ(K,J) = RZ(J,K)
            SSPZ(K,J) = SSPZ(J,K)
  660    CONTINUE
  680 CONTINUE
      NCASES = INT(XCASES+0.5D0)
      IFAIL = 0
      RETURN
  700 IERROR = 4
      GO TO 760
  720 NCASES = 0
      IERROR = 6
      GO TO 760
  740 NCASES = 1
      IERROR = 7
  760 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
