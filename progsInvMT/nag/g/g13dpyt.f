      SUBROUTINE G13DPY(MEAN,WEIGHT,M,ISX,Q,LDQ,IP,X,NY,Y,WT,RSS,WK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     ADDS (UPDATE='A') AN OBSERVATION FROM A MULTIVARIATE LINEAR
C     REGRESSION RETURNING THE UPDATED RSS AND NEW R MATRIX AND Q'Y
C     MATRIX STORED IN Q.
C     BASED ON G02DCF.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  WT
      INTEGER           IP, LDQ, M, NY
      CHARACTER         MEAN, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,IP+NY), RSS(NY,NY), WK(3*IP+NY), X(*),
     *                  Y(NY)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  QYI, SQRWT
      INTEGER           I, J, K
C     .. External Subroutines ..
      EXTERNAL          F06QQF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
         K = 1
      ELSE
         K = 0
      END IF
      DO 20 I = 1, M
         IF (ISX(I).GT.0) THEN
            K = K + 1
         END IF
   20 CONTINUE
      IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
         SQRWT = SQRT(WT)
         DO 40 I = 1, NY
            WK(3*IP+I) = SQRWT*Y(I)
   40    CONTINUE
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            K = 1
            WK(2*IP+1) = SQRWT
         ELSE
            K = 0
         END IF
         DO 60 I = 1, M
            IF (ISX(I).GT.0) THEN
               K = K + 1
               WK(2*IP+K) = SQRWT*X(I)
            END IF
   60    CONTINUE
      ELSE
         DO 80 I = 1, NY
            WK(3*IP+I) = Y(I)
   80    CONTINUE
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            K = 1
            WK(2*IP+1) = 1.0D0
         ELSE
            K = 0
         END IF
         DO 100 I = 1, M
            IF (ISX(I).GT.0) THEN
               K = K + 1
               WK(2*IP+K) = X(I)
            END IF
  100    CONTINUE
      END IF
C
C           UPDATE R
C
      CALL F06QQF(IP,1.0D0,WK(2*IP+1),1,Q(1,NY+1),LDQ,WK(1),WK(IP+1))
C
C           UPDATE Q'Y
C
      DO 140 I = 1, IP
         DO 120 J = 1, NY
            QYI = Q(I,J)
            Q(I,J) = WK(IP+I)*WK(3*IP+J) + WK(I)*QYI
            WK(3*IP+J) = WK(I)*WK(3*IP+J) - WK(IP+I)*QYI
  120    CONTINUE
  140 CONTINUE
C
C           UPDATE RSS
C
      DO 180 I = 1, NY
         DO 160 J = 1, I
            RSS(J,I) = WK(3*IP+J)*WK(3*IP+I) + RSS(J,I)
  160    CONTINUE
  180 CONTINUE
      RETURN
      END
