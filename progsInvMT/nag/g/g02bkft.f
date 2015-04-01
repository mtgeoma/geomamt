      SUBROUTINE G02BKF(N,M,X,IX,NVARS,KVAR,XBAR,STD,SSPZ,ISSPZ,RZ,IRZ,
     *                  IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 10 REVISED. IER-381 (JUN 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BKF
C     WRITTEN 17. 8.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES MEANS AND STANDARD DEVIATIONS OF VARIABLES, SUMS OF
C     SQUARES AND CROSS-PRODUCTS ABOUT ZERO AND CORRELATION-LIKE
C     COEFFICIENTS FOR A SET OF DATA IN SPECIFIED COLUNNS OF THE
C     ARRAY X
C
C     USES NAG ERROR ROUTINE P01AAF
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE-DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH MAY BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BKF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IRZ, ISSPZ, IX, M, N, NVARS
C     .. Array Arguments ..
      DOUBLE PRECISION  RZ(IRZ,NVARS), SSPZ(ISSPZ,NVARS), STD(NVARS),
     *                  X(IX,M), XBAR(NVARS)
      INTEGER           KVAR(NVARS)
C     .. Local Scalars ..
      DOUBLE PRECISION  FN, S, XM, XN
      INTEGER           I, IERROR, J, JM, JV, K, KV
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      IERROR = 0
      IF (IX.LT.N .OR. ISSPZ.LT.NVARS .OR. IRZ.LT.NVARS) IERROR = 3
      IF (NVARS.LT.2 .OR. NVARS.GT.M) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 300
   20 DO 40 I = 1, NVARS
         IF (KVAR(I).LT.1 .OR. KVAR(I).GT.M) GO TO 280
   40 CONTINUE
      FN = DBLE(N)
      DO 260 J = 1, NVARS
         JV = KVAR(J)
         S = 0.0D0
         DO 60 I = 1, N
            S = S + X(I,JV)
   60    CONTINUE
         XM = S/FN
         XBAR(J) = XM
         S = 0.0D0
         DO 80 I = 1, N
            XN = X(I,JV) - XM
            S = S + XN*XN
   80    CONTINUE
         SSPZ(J,J) = S + XM*XM*FN
         IF (S) 100, 100, 120
  100    STD(J) = 0.0D0
         RZ(J,J) = 0.0D0
         IF (SSPZ(J,J).GT.0.0D0) RZ(J,J) = 1.0D0
         GO TO 140
  120    STD(J) = SQRT(S/(FN-1.0D0))
         RZ(J,J) = 1.0D0
  140    IF (J.EQ.1) GO TO 260
         JM = J - 1
         DO 240 K = 1, JM
            S = 0.0D0
            KV = KVAR(K)
            DO 160 I = 1, N
               S = S + X(I,JV)*X(I,KV)
  160       CONTINUE
            SSPZ(J,K) = S
            SSPZ(K,J) = S
            S = SSPZ(J,J)*SSPZ(K,K)
            IF (S) 180, 180, 200
  180       RZ(J,K) = 0.0D0
            GO TO 220
  200       RZ(J,K) = SSPZ(J,K)/SQRT(S)
  220       RZ(K,J) = RZ(J,K)
  240    CONTINUE
  260 CONTINUE
      IFAIL = 0
      RETURN
  280 IERROR = 4
  300 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
