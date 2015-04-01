      SUBROUTINE G10ABZ(WEIGHT,X,AVH,Y,WT,AVDY,N,WWT,U,C,LDC,R,T,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Initializes the arrays C, R and T for one dimensional cubic
C     Smoothing spline fitting by subroutine G10ABY.  The values
C     WT(I) are scaled so that the sum of their squares is N
C     and the average of the differences X(I+1) - X(I) is calculated
C     in AVH in order to avoid underflow and overflow problems in
C     G10ABY.
C
C     Subroutine sets IERROR if elements of X are non-increasing
C     or if WT(I) is not positive for some I.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AVDY, AVH
      INTEGER           IERROR, LDC, N
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,3), R(N+2,3), T(N+2,2), U(N), WT(*),
     *                  WWT(N), X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, E, G, H, ZERO
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
C
C     Get average X spacing in AVH
C
      G = ZERO
      DO 20 I = 1, N - 1
         H = X(I+1) - X(I)
         IF (H.LE.ZERO) THEN
            IERROR = 3
            GO TO 140
         END IF
         G = G + H
   20 CONTINUE
      AVH = G/(N-1)
C
C     Scale relative weights
C
      IF ((WEIGHT.EQ.'W') .OR. (WEIGHT.EQ.'w')) THEN
         G = ZERO
         DO 40 I = 1, N
            IF (WT(I).LE.ZERO) THEN
               IERROR = 2
               GO TO 140
            END IF
            G = G + 1.0D0/WT(I)
   40    CONTINUE
         AVDY = SQRT(G/N)
         DO 60 I = 1, N
            WWT(I) = 1.0D0/(SQRT(WT(I))*AVDY)
   60    CONTINUE
      ELSE
         AVDY = 1.0D0
         DO 80 I = 1, N
            WWT(I) = 1.0D0
   80    CONTINUE
      END IF
C
C     Initialize H,F
C
      H = (X(2)-X(1))/AVH
      E = (Y(2)-Y(1))/H
C
C     Calculate T,R
C
      DO 100 I = 2, N - 1
         G = H
         H = (X(I+1)-X(I))/AVH
         D = E
         E = (Y(I+1)-Y(I))/H
         U(I) = E - D
         T(I+1,1) = 2.0D0*(G+H)/3.0D0
         T(I+1,2) = H/3.0D0
         R(I+1,3) = WWT(I-1)/G
         R(I+1,1) = WWT(I+1)/H
         R(I+1,2) = -WWT(I)/G - WWT(I)/H
  100 CONTINUE
C
C     Calculate C = R'*R
C
      R(N+1,2) = ZERO
      R(N+1,3) = ZERO
      R(N+2,3) = ZERO
      DO 120 I = 2, N - 1
         C(I,1) = R(I+1,1)*R(I+1,1) + R(I+1,2)*R(I+1,2) + R(I+1,3)
     *            *R(I+1,3)
         C(I,2) = R(I+1,1)*R(I+2,2) + R(I+1,2)*R(I+2,3)
         C(I,3) = R(I+1,1)*R(I+3,3)
  120 CONTINUE
  140 CONTINUE
      RETURN
      END
