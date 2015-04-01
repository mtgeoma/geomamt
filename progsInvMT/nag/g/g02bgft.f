      SUBROUTINE G02BGF(N,M,X,IX,NVARS,KVAR,XBAR,STD,SSP,ISSP,R,IR,
     *                  IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BGF
C     WRITTEN 17. 7.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES MEANS, STANDARD DEVIATIONS, SUMS OF SQUARES AND
C     CROSS-PRODUCTS OF DEVIATIONS FROM MEANS, AND PEARSON PRODUCT-
C     MOMENT CORRELATION COEFFICIENTS FOR A SET OF DATA IN
C     SPECIFIED COLUNNS OF THE ARRAY X
C
C     USES NAG ERROR ROUTINE P01AAF
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH MAY BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BGF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IR, ISSP, IX, M, N, NVARS
C     .. Array Arguments ..
      DOUBLE PRECISION  R(IR,NVARS), SSP(ISSP,NVARS), STD(NVARS),
     *                  X(IX,M), XBAR(NVARS)
      INTEGER           KVAR(NVARS)
C     .. Local Scalars ..
      DOUBLE PRECISION  FN, S
      INTEGER           I, IERROR, J, JV, K, KV
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      IERROR = 0
      IF (IX.LT.N .OR. ISSP.LT.NVARS .OR. IR.LT.NVARS) IERROR = 3
      IF (NVARS.LT.2 .OR. NVARS.GT.M) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 340
   20 DO 40 I = 1, NVARS
         IF (KVAR(I).LT.1 .OR. KVAR(I).GT.M) GO TO 320
   40 CONTINUE
      FN = DBLE(N)
      DO 120 J = 1, NVARS
         JV = KVAR(J)
         S = 0.0D0
         DO 60 I = 1, N
            S = S + X(I,JV)
   60    CONTINUE
         XBAR(J) = S/FN
         DO 100 K = 1, J
            KV = KVAR(K)
            S = 0.0D0
            DO 80 I = 1, N
               S = S + (X(I,JV)-XBAR(J))*(X(I,KV)-XBAR(K))
   80       CONTINUE
            SSP(J,K) = S
            SSP(K,J) = S
  100    CONTINUE
  120 CONTINUE
      DO 280 J = 1, NVARS
         S = SSP(J,J)
         IF (S) 140, 140, 160
  140    STD(J) = 0.0D0
         GO TO 180
  160    STD(J) = SQRT(S)
  180    DO 260 K = 1, J
            S = STD(J)*STD(K)
            IF (S) 200, 200, 220
  200       R(J,K) = 0.0D0
            GO TO 240
  220       R(J,K) = SSP(J,K)/S
  240       R(K,J) = R(J,K)
  260    CONTINUE
  280 CONTINUE
      S = SQRT(FN-1.0D0)
      DO 300 J = 1, NVARS
         STD(J) = STD(J)/S
  300 CONTINUE
      IFAIL = 0
      RETURN
  320 IERROR = 4
  340 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
