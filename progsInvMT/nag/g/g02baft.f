      SUBROUTINE G02BAF(N,M,X,IX,XBAR,STD,SSP,ISSP,R,IR,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BAF
C     WRITTEN 16. 7.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C     (REVAMP OF G02AAF ALLOWING MORE FLEXIBLE ARRAY DIMENSIONS)
C
C     COMPUTES MEANS, STANDARD DEVIATIONS, SUMS OF SQUARES AND
C     CROSS-PRODUCTS OF DEVIATIONS FROM MEANS, AND PEARSON PRODUCT-
C     MOMENT CORRELATION COEFFICIENTS FOR A SET OF DATA IN THE
C     ARRAY X.
C
C     USES NAG ERROR ROUTINE P01AAF
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IR, ISSP, IX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  R(IR,M), SSP(ISSP,M), STD(M), X(IX,M), XBAR(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  FN, S
      INTEGER           I, IERROR, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      IERROR = 0
      IF (IX.LT.N .OR. ISSP.LT.M .OR. IR.LT.M) IERROR = 3
      IF (M.LT.2) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 300
   20 FN = DBLE(N)
      DO 100 J = 1, M
         S = 0.0D0
         DO 40 I = 1, N
            S = S + X(I,J)
   40    CONTINUE
         XBAR(J) = S/FN
         DO 80 K = 1, J
            S = 0.0D0
            DO 60 I = 1, N
               S = S + (X(I,J)-XBAR(J))*(X(I,K)-XBAR(K))
   60       CONTINUE
            SSP(J,K) = S
            SSP(K,J) = S
   80    CONTINUE
  100 CONTINUE
      DO 260 J = 1, M
         S = SSP(J,J)
         IF (S) 120, 120, 140
  120    STD(J) = 0.0D0
         GO TO 160
  140    STD(J) = SQRT(S)
  160    DO 240 K = 1, J
            S = STD(J)*STD(K)
            IF (S) 180, 180, 200
  180       R(J,K) = 0.0D0
            GO TO 220
  200       R(J,K) = SSP(J,K)/S
  220       R(K,J) = R(J,K)
  240    CONTINUE
  260 CONTINUE
      S = SQRT(FN-1.0D0)
      DO 280 J = 1, M
         STD(J) = STD(J)/S
  280 CONTINUE
      IFAIL = 0
      RETURN
  300 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
