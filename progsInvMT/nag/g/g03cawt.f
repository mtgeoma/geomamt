      SUBROUTINE G03CAW(WEIGHT,N,X,LDX,M,ISX,IVAR,WT,T,V,LDV,S,E)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes standardized data values
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IVAR, LDV, LDX, M, N
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  E(IVAR), S(IVAR), V(LDV,*), WT(*), X(LDX,M)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  SCALE, SUM, XBAR
      INTEGER           I, J, K
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
C     Compute means
C
      K = 0
      DO 60 J = 1, M
         IF (ISX(J).GT.0) THEN
            K = K + 1
            SUM = 0.0D0
            IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
               DO 20 I = 1, N
                  SUM = SUM + X(I,J)
   20          CONTINUE
               E(K) = SUM/T
            ELSE
               DO 40 I = 1, N
                  IF (WT(I).GT.0.0D0) SUM = SUM + WT(I)*X(I,J)
   40          CONTINUE
               E(K) = SUM/T
            END IF
         END IF
   60 CONTINUE
C
C     Compute variances and standardize
C
      K = 0
      IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
C
C        Weighted case
C
         DO 140 J = 1, M
            IF (ISX(J).GT.0) THEN
               K = K + 1
               XBAR = E(K)
               SUM = 0.0D0
               DO 80 I = 1, N
                  SUM = SUM + (X(I,J)-XBAR)**2
   80          CONTINUE
               IF (SUM.GT.0.0D0) THEN
                  SCALE = 1.0D0/SQRT(SUM)
                  DO 100 I = 1, N
                     V(I,K) = (X(I,J)-XBAR)*SCALE
  100             CONTINUE
               ELSE
                  DO 120 I = 1, N
                     V(I,K) = 0.0D0
  120             CONTINUE
               END IF
               S(K) = SUM
            END IF
  140    CONTINUE
      ELSE
C
C        unweighted case
C
         DO 220 J = 1, M
            IF (ISX(J).GT.0) THEN
               K = K + 1
               XBAR = E(K)
               SUM = 0.0D0
               DO 160 I = 1, N
                  IF (WT(I).GT.0.0D0) SUM = SUM + WT(I)*((X(I,J)-XBAR)
     *                                      **2)
  160          CONTINUE
               IF (SUM.GT.0.0D0) THEN
                  SCALE = 1.0D0/SQRT(SUM)
                  DO 180 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        V(I,K) = V(I,IVAR)*(X(I,J)-XBAR)*SCALE
                     ELSE
                        V(I,K) = 0.0D0
                     END IF
  180             CONTINUE
               ELSE
                  DO 200 I = 1, N
                     V(I,K) = 0.0D0
  200             CONTINUE
               END IF
               S(K) = SUM
            END IF
  220    CONTINUE
      END IF
      RETURN
      END
