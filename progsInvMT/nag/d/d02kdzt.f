      SUBROUTINE D02KDZ(X,H,N,Y,FCN,W,IW1,IW2,COEFFN,COEFF1,M,ARR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-227 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     COEFF1, COEFFN, FCN
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, X
      INTEGER           IW1, IW2, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ARR(M), W(IW1,IW2), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          COEFF1, COEFFN, FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, S
      INTEGER           I, N1
C     .. Local Arrays ..
      DOUBLE PRECISION  C(10)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      EPS = X02AJF()
      C(1) = 1.D0/6.D0
      C(2) = 1.D0/3.D0
      C(3) = 0.125D0
      C(4) = 0.375D0
      C(5) = 0.5D0
      C(6) = 1.5D0
      C(7) = 2.0D0
      C(8) = 2.D0/3.D0
      C(9) = 0.2D0
      C(10) = 4.D0/3.D0
      N1 = 6
      IF (IW2.EQ.4) N1 = 4
      IF (IW2.EQ.6) N1 = 5
      DO 20 I = 1, N
         W(I,3) = Y(I) + C(2)*H*W(I,1)
   20 CONTINUE
      CALL FCN(N,X+C(2)*H,W(1,3),W(1,N1),COEFFN,COEFF1,M,ARR)
      DO 40 I = 1, N
         IF (N1.EQ.5) W(I,4) = 0.5D0*(W(I,N1)-W(I,1))*C(2)*H
         W(I,3) = Y(I) + C(1)*H*(W(I,1)+W(I,N1))
   40 CONTINUE
      CALL FCN(N,X+C(2)*H,W(1,3),W(1,N1),COEFFN,COEFF1,M,ARR)
      DO 60 I = 1, N
         W(I,2) = Y(I) + H*(C(3)*W(I,1)+C(4)*W(I,N1))
   60 CONTINUE
      CALL FCN(N,X+C(5)*H,W(1,2),W(1,3),COEFFN,COEFF1,M,ARR)
      DO 100 I = 1, N
         IF (N1.EQ.4) GO TO 80
         W(I,IW2) = -C(2)*W(I,1) - C(10)*W(I,3) + C(6)*W(I,N1)
   80    W(I,2) = Y(I) + H*(C(5)*W(I,1)-C(6)*W(I,N1)+C(7)*W(I,3))
  100 CONTINUE
      CALL FCN(N,X+H,W(1,2),W(1,N1),COEFFN,COEFF1,M,ARR)
      DO 140 I = 1, N
         W(I,2) = Y(I)
         Y(I) = Y(I) + H*(C(1)*(W(I,1)+W(I,N1))+C(8)*W(I,3))
         IF (N1.EQ.4) GO TO 120
         W(I,3) = W(I,N1)
         W(I,N1) = C(9)*H*(W(I,IW2)+C(1)*W(I,N1))
         S = 0.D0
         IF (ABS(W(I,N1)).LE.30.D0*C(9)*EPS*ABS(H)*MAX(ABS(W(I,IW2))
     *       ,C(1)*ABS(W(I,3)))) S = 1.D0
         W(I,IW2) = S
  120    W(I,3) = W(I,1)
  140 CONTINUE
      RETURN
      END
