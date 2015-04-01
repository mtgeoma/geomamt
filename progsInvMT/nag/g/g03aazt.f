      SUBROUTINE G03AAZ(MATRIX,WEIGHT,N,X,LDX,M,ISX,IVAR,WT,T,V,LDV,S,E)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IVAR, LDV, LDX, M, N
      CHARACTER         MATRIX, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  E(IVAR), S(*), V(LDV,*), WT(*), X(LDX,M)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  SCALE, SUM, XBAR
      INTEGER           I, J, K
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      K = 0
      DO 60 J = 1, M
         IF (((MATRIX.NE.'N' .AND. MATRIX.NE.'n') .AND. (ISX(J).GT.0))
     *        .OR. ((MATRIX.EQ.'N' .OR. MATRIX.EQ.'n') .AND. (ISX(J)
     *       .LT.0))) THEN
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
      K = 0
      IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
         IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') THEN
            DO 120 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  SUM = 0.0D0
                  DO 80 I = 1, N
                     SUM = SUM + (X(I,J)-XBAR)**2
   80             CONTINUE
                  IF (SUM.GT.0.0D0) THEN
                     SUM = SQRT(SUM)
                     DO 100 I = 1, N
                        V(I,K) = (X(I,J)-XBAR)/SUM
  100                CONTINUE
                  END IF
                  S(J) = SUM
               END IF
  120       CONTINUE
         ELSE IF (MATRIX.EQ.'S' .OR. MATRIX.EQ.'s') THEN
            DO 160 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  SCALE = S(J)
                  DO 140 I = 1, N
                     V(I,K) = (X(I,J)-XBAR)/SCALE
  140             CONTINUE
               END IF
  160       CONTINUE
         ELSE IF (MATRIX.EQ.'U' .OR. MATRIX.EQ.'u' .OR. MATRIX.EQ.
     *            'V' .OR. MATRIX.EQ.'v') THEN
            DO 200 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  DO 180 I = 1, N
                     V(I,K) = X(I,J) - XBAR
  180             CONTINUE
               END IF
  200       CONTINUE
         ELSE IF (MATRIX.EQ.'N' .OR. MATRIX.EQ.'n') THEN
            DO 240 J = 1, M
               IF (ISX(J).LT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  DO 220 I = 1, N
                     V(I,K) = X(I,J) - XBAR
  220             CONTINUE
               END IF
  240       CONTINUE
         END IF
      ELSE
         IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') THEN
            DO 300 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  SUM = 0.0D0
                  DO 260 I = 1, N
                     IF (WT(I).GT.0.0D0) SUM = SUM + WT(I)*((X(I,J)
     *                   -XBAR)**2)
  260             CONTINUE
                  IF (SUM.GT.0.0D0) THEN
                     SUM = SQRT(SUM)
                     DO 280 I = 1, N
                        IF (WT(I).GT.0.0D0) THEN
                           V(I,K) = V(I,IVAR)*(X(I,J)-XBAR)/SUM
                        ELSE
                           V(I,K) = 0.0D0
                        END IF
  280                CONTINUE
                  END IF
                  S(J) = SUM
               END IF
  300       CONTINUE
         ELSE IF (MATRIX.EQ.'S' .OR. MATRIX.EQ.'s') THEN
            DO 340 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  SCALE = S(J)
                  DO 320 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        V(I,K) = V(I,IVAR)*(X(I,J)-XBAR)/SCALE
                     ELSE
                        V(I,K) = 0.0D0
                     END IF
  320             CONTINUE
               END IF
  340       CONTINUE
         ELSE IF (MATRIX.EQ.'U' .OR. MATRIX.EQ.'u' .OR. MATRIX.EQ.
     *            'V' .OR. MATRIX.EQ.'v') THEN
            DO 380 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  DO 360 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        V(I,K) = V(I,IVAR)*(X(I,J)-XBAR)
                     ELSE
                        V(I,K) = 0.0D0
                     END IF
  360             CONTINUE
               END IF
  380       CONTINUE
         ELSE IF (MATRIX.EQ.'N' .OR. MATRIX.EQ.'n') THEN
            DO 420 J = 1, M
               IF (ISX(J).LT.0) THEN
                  K = K + 1
                  XBAR = E(K)
                  DO 400 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        V(I,K) = V(I,IVAR)*(X(I,J)-XBAR)
                     ELSE
                        V(I,K) = 0.0D0
                     END IF
  400             CONTINUE
               END IF
  420       CONTINUE
         END IF
      END IF
      RETURN
      END
