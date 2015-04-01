      SUBROUTINE G03FCY(K,X,XHAT,TOL,W)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Average elements of X as necessary to achieve monotonicity.
C     Uses the up-and-down blocks algorithm s given by Cran (Applied
C     Statistics, 1980) but modified to avoid array shifts. Algorithm
C     only works for positive X.
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           K
C     .. Array Arguments ..
      DOUBLE PRECISION  W(K), X(K), XHAT(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  S, WW
      INTEGER           I, I1, IM1, INC, J
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      DO 20 I = 1, K
         W(I) = 1.0D0
   20 CONTINUE
      I = 1
   40 CONTINUE
      IF (I.LT.K) THEN
         INC = 1
   60    CONTINUE
         IF (X(I+INC).LT.0.0D0) THEN
            IF (I+INC.EQ.K) GO TO 80
            INC = INC + 1
            GO TO 60
         END IF
         IF (X(I).GT.X(I+INC)) THEN
            GO TO 180
         ELSE IF (I.EQ.1) THEN
            GO TO 260
         END IF
      END IF
   80 CONTINUE
      INC = 1
  100 CONTINUE
      IF (X(I-INC).LT.0.0D0) THEN
         INC = INC + 1
         IF (I-INC.LT.1) THEN
            INC = 1
  120       CONTINUE
            IF (X(I+INC).LT.0.0D0) THEN
               IF (I+INC.EQ.K) GO TO 280
               INC = INC + 1
               GO TO 120
            END IF
            GO TO 260
         END IF
         GO TO 100
      END IF
      IF (X(I-INC).LE.X(I)) THEN
         IF (I.LT.K) THEN
            INC = 1
  140       CONTINUE
            IF (X(I+INC).LT.0.0D0) THEN
               IF (I+INC.EQ.K) GO TO 280
               INC = INC + 1
               GO TO 140
            END IF
            GO TO 260
         ELSE
            GO TO 280
         END IF
      END IF
C
C     Pool the active block with the next lower block
C
  160 IM1 = I - INC
      WW = W(IM1) + W(I)
      X(IM1) = (W(IM1)*X(IM1)+W(I)*X(I))/WW
      W(IM1) = WW
      X(I) = -1.0D0
      I = IM1
      GO TO 40
C
C     Pool the active block with the next higher block
C
  180 I1 = I + INC
      WW = W(I) + W(I1)
      X(I) = (W(I)*X(I)+W(I1)*X(I1))/WW
      W(I) = WW
      X(I1) = -1.0D0
      IF (I.GT.1) THEN
         INC = 1
  200    CONTINUE
         IF (X(I-INC).LT.0.0D0) THEN
            IF (I-INC.EQ.1) GO TO 220
            INC = INC + 1
            GO TO 200
         END IF
         IF (X(I-INC).GT.X(I)) GO TO 160
      END IF
  220 CONTINUE
      IF (I.GE.K) GO TO 280
      INC = 1
  240 CONTINUE
      IF (X(I+INC).LT.0.0D0) THEN
         IF (I+INC.EQ.K) GO TO 280
         INC = INC + 1
         GO TO 240
      END IF
      IF (X(I).GT.X(I+INC)) GO TO 180
  260 I = I + INC
      GO TO 220
C
C     Obtain the amalgamated means XHAT from the working array X
C
  280 CONTINUE
      I1 = 1
      DO 340 I = 1, K
         IF (X(I).GT.0.0D0) THEN
            S = 0.0D0
            DO 300 J = I1, K
               S = S + 1.0D0
               XHAT(J) = X(I)
               IF (ABS(S-W(I)).LT.TOL) GO TO 320
  300       CONTINUE
  320       I1 = J + 1
         END IF
  340 CONTINUE
      RETURN
      END
