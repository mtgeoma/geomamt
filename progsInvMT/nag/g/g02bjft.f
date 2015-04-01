      SUBROUTINE G02BJF(N,M,X,IX,MISS,XMISS,NVARS,KVAR,XBAR,STD,SSP,
     *                  ISSP,R,IR,NCASES,COUNT,IC,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BJF
C     WRITTEN 17. 8.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES MEANS, STANDARD DEVIATIONS, SUMS OF SQUARES AND
C     CROSS-PRODUCTS OF DEVIATIONS FROM MEANS, AND PEARSON PRODUCT-
C     MOMENT CORRELATION COEFFICIENTS FOR A SET OF DATA IN
C     SPECIFIED COLUMNS OF THE ARRAY X, OMITTING CASES WITH MISSING
C     VALUES FROM ONLY THOSE CALCULATIONS INVOLVING THE VARIABLES
C     FOR WHICH THE VALUES ARE MISSING.
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
      PARAMETER         (SRNAME='G02BJF')
C     .. Scalar Arguments ..
      INTEGER           IC, IFAIL, IR, ISSP, IX, M, N, NCASES, NVARS
C     .. Array Arguments ..
      DOUBLE PRECISION  COUNT(IC,NVARS), R(IR,NVARS), SSP(ISSP,NVARS),
     *                  STD(NVARS), X(IX,M), XBAR(NVARS), XMISS(M)
      INTEGER           KVAR(NVARS), MISS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, XJ, XK, XN, XNM
      INTEGER           I, IERROR, J, JP, JV, K, KV
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
      IF (IX.LT.N .OR. ISSP.LT.NVARS .OR. IR.LT.NVARS .OR. IC.LT.NVARS)
     *    IERROR = 3
      IF (NVARS.LT.2 .OR. NVARS.GT.M) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 540
   20 DO 40 I = 1, NVARS
         IF (KVAR(I).LT.1 .OR. KVAR(I).GT.M) GO TO 500
   40 CONTINUE
C
C
      DO 80 J = 1, NVARS
         DO 60 K = 1, NVARS
            SSP(K,J) = 0.0D0
            COUNT(K,J) = 0.0D0
   60    CONTINUE
         R(J,J) = 0.0D0
   80 CONTINUE
C
C     ITERATIVE PROCESS THROUGH CASES BUILDING UP SUMS OF SQUARES
C     AND CROSS-PRODUCTS OF DEVIATIONS
C
      DO 320 I = 1, N
         DO 300 J = 1, NVARS
            JV = KVAR(J)
            IF (MISS(JV)) 120, 120, 100
  100       IF (ABS(X(I,JV)-XMISS(JV)).LE.ABS(ACC*XMISS(JV)))
     *          GO TO 300
  120       IF (COUNT(J,J)) 140, 140, 160
  140       COUNT(J,J) = 1.0D0
            R(J,J) = X(I,JV)
            GO TO 180
  160       XNM = COUNT(J,J)
            XN = XNM + 1.0D0
            XJ = (X(I,JV)-R(J,J))/XN
            SSP(J,J) = SSP(J,J) + XJ*XJ*XN*XNM
            R(J,J) = R(J,J) + XJ
            COUNT(J,J) = XN
  180       IF (J.EQ.NVARS) GO TO 300
            JP = J + 1
            DO 280 K = JP, NVARS
               KV = KVAR(K)
               IF (MISS(KV)) 220, 220, 200
  200          IF (ABS(X(I,KV)-XMISS(KV)).LE.ABS(ACC*XMISS(KV)))
     *             GO TO 280
  220          IF (COUNT(K,J)) 240, 240, 260
  240          COUNT(K,J) = 1.0D0
               R(K,J) = X(I,KV)
               R(J,K) = X(I,JV)
               GO TO 280
  260          XNM = COUNT(K,J)
               XN = XNM + 1.0D0
               XJ = (X(I,JV)-R(J,K))/XN
               XK = (X(I,KV)-R(K,J))/XN
               COUNT(K,J) = XN
               COUNT(J,K) = COUNT(J,K) + XJ*XJ*XNM*XN
               SSP(K,J) = SSP(K,J) + XK*XK*XNM*XN
               SSP(J,K) = SSP(J,K) + XJ*XK*XNM*XN
               R(K,J) = R(K,J) + XK
               R(J,K) = R(J,K) + XJ
  280       CONTINUE
  300    CONTINUE
  320 CONTINUE
C
C     TIDY UP MATRICES TO PROVIDE PROMISED OUTPUT
C
      DO 440 J = 1, NVARS
         XBAR(J) = R(J,J)
         XJ = SSP(J,J)
         XNM = COUNT(J,J) - 1.0D0
         IF (XJ.LE.0.0D0 .OR. XNM.LT.0.5D0) GO TO 340
         STD(J) = SQRT(XJ/XNM)
         R(J,J) = 1.0D0
         GO TO 360
  340    STD(J) = 0.0D0
         R(J,J) = 0.0D0
  360    IF (J.EQ.NVARS) GO TO 440
         JP = J + 1
         DO 420 K = JP, NVARS
            XJ = SSP(K,J)*COUNT(J,K)
            IF (XJ.LE.0.0D0 .OR. COUNT(K,J).LT.1.5D0) GO TO 380
            R(J,K) = SSP(J,K)/SQRT(XJ)
            GO TO 400
  380       R(J,K) = 0.0D0
  400       R(K,J) = R(J,K)
            SSP(K,J) = SSP(J,K)
            COUNT(J,K) = COUNT(K,J)
  420    CONTINUE
  440 CONTINUE
C
C     CHECK ON NUMBERS OF CASES
C
      XN = DBLE(N)
      DO 480 J = 1, NVARS
         DO 460 K = J, NVARS
            IF (COUNT(K,J).LT.XN) XN = COUNT(K,J)
  460    CONTINUE
  480 CONTINUE
      NCASES = INT(XN+0.5D0)
      IF (NCASES.LT.2) GO TO 520
      IFAIL = 0
      RETURN
  500 IERROR = 4
      GO TO 540
  520 IERROR = 5
  540 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
