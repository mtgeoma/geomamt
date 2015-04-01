      SUBROUTINE F02BFF(C,B,BETA,N,M1,M2,MM12,EPS1,RELFEH,EPS2,IZ,X,WU)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-515 (AUG 1986).
C     WRITTEN BY W.PHILLIPS     1ST OCTOBER 1975
C     OXFORD UNIVERSITY COMPUTING LABORATORY.
C     THIS ROUTINE REPLACES F02AZF.
C
C     BISECT
C     C IS THE DIAGONAL, B THE SUB-DIAGONAL AND BETA THE SQUARED
C     SUBDIAGONAL OF A SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N. THE
C     EIGENVALUES LAMBDA(M1)......LAMBDA(M2), WHERE M2 IS NOT
C     LESS THAN M1 AND LAMBDA(I+1) IS NOT LESS THAN LAMBDA(I),
C     ARE CALCULATED BY THE METHOD OF BISECTION AND STORED IN THE
C     VECTOR X. BISECTION IS CONTINUED UNTIL THE UPPER AND LOWER
C     BOUNDS FOR AN EIGENVALUE DIFFER BY LESS THAN EPS1, UNLESS AT
C     SOME EARLIER STAGE, THE UPPER AND LOWER BOUNDS DIFFER ONLY
C     IN THE LEAST SIGNIFICANT DIGITS. EPS2 GIVES AN EXTREME UPPER
C     BOUND FOR THE ERROR IN ANY EIGENVALUE, BUT FOR CERTAIN TYPES
C     OF MATRICES THE SMALL EIGENVALUES ARE DETERMINED TO A VERY
C     MUCH HIGHER ACCURACY. IN THIS CASE, EPS1 SHOULD BE SET EQUAL
C     TO THE ERROR TO BE TOLERATED IN THE SMALLEST EIGENVALUE.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS1, EPS2, RELFEH
      INTEGER           IZ, M1, M2, MM12, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), BETA(N), C(N), WU(MM12), X(MM12)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HALF, ONE, Q, SEVEN, TWO, X0, X1, XMAX, XMIN,
     *                  XU, Y, ZERO
      INTEGER           I, IA, IAM, II, IM, K, KK, KM, N1
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ONE/1.0D0/, ZERO/0.0D0/, TWO/2.0D0/,
     *                  HALF/0.5D0/, SEVEN/7.0D0/
C     .. Executable Statements ..
C     CALCULATION OF XMIN, XMAX
      BETA(1) = ZERO
      B(1) = ZERO
      XMIN = C(N) - ABS(B(N))
      XMAX = C(N) + ABS(B(N))
      I = N
      N1 = N - 1
      IF (N1.LT.1) GO TO 40
      DO 20 II = 1, N1
         I = I - 1
         H = ABS(B(I)) + ABS(B(I+1))
         IF ((C(I)+H).GT.XMAX) XMAX = C(I) + H
         IF ((C(I)-H).LT.XMIN) XMIN = C(I) - H
   20 CONTINUE
   40 Y = XMAX
      IF ((XMIN+XMAX).LE.ZERO) Y = -XMIN
      EPS2 = RELFEH*Y
      IF (EPS1.LE.ZERO) EPS1 = EPS2
      EPS2 = HALF*EPS1 + SEVEN*EPS2
C     INNER BLOCK
      X0 = XMAX
      IF (M2.LT.M1) GO TO 80
      DO 60 I = M1, M2
         II = I - M1 + 1
         X(II) = XMAX
         WU(II) = XMIN
   60 CONTINUE
   80 IZ = 0
C     LOOP FOR K-TH EIGENVALUE
      IF (M2.LT.M1) RETURN
      K = M2 + 1
      DO 240 KK = M1, M2
         K = K - 1
         XU = XMIN
         I = K + 1
         IF (K.LT.M1) GO TO 120
         DO 100 II = M1, K
            I = I - 1
            IM = I - M1 + 1
            IF (XU.GE.WU(IM)) GO TO 100
            XU = WU(IM)
            GO TO 120
  100    CONTINUE
  120    KM = K - M1 + 1
         IF (X0.GT.X(KM)) X0 = X(KM)
  140    X1 = (XU+X0)*HALF
         IF ((X0-XU).LE.(TWO*RELFEH*(ABS(XU)+ABS(X0))+EPS1)) GO TO 220
         IZ = IZ + 1
C        STURMS SEQUENCE
         IA = 0
         Q = ONE
         DO 160 I = 1, N
            IF (Q.EQ.ZERO) Y = ABS(B(I)/RELFEH)
            IF (Q.NE.ZERO) Y = BETA(I)/Q
            Q = C(I) - X1 - Y
            IF (Q.LT.ZERO) IA = IA + 1
  160    CONTINUE
         IF (IA.GE.K) GO TO 200
         IF (IA.GE.M1) GO TO 180
         XU = X1
         WU(1) = X1
         GO TO 140
  180    IAM = IA - M1 + 2
         XU = X1
         WU(IAM) = X1
         IAM = IA - M1 + 1
         IF (X(IAM).GT.X1) X(IAM) = X1
         GO TO 140
  200    X0 = X1
         GO TO 140
  220    KM = K - M1 + 1
         X(KM) = (X0+XU)*HALF
  240 CONTINUE
      RETURN
      END
