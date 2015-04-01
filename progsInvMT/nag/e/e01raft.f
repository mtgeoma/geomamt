      SUBROUTINE E01RAF(N,X,F,M,A,U,IW,IFAIL)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     E01RAF PRODUCES, FROM A SET OF FUNCTION VALUES
C     AND CORRESPONDING ABSCISSAE, THE COEFFICIENTS OF
C     AN INTERPOLATING RATIONAL FUNCTION EXPRESSED IN
C     CONTINUED FRACTION FORM
C
C     ON ENTRY
C     N  -INTEGER-  NUMBER OF INTERPOLATION POINTS
C     X  -REAL ARRAY-  ABSCISSA POINTS
C     F  -REAL ARRAY-  FUNCTION VALUES
C     ON EXIT
C     M  -INTEGER-  NUMBER OF THIELE COEFFICIENTS
C     U  -REAL ARRAY-  REORDERED ABSCISSA POINTS
C     A  -REAL ARRAY-  THIELE COEFFICIENTS
C     IW  -INTEGER ARRAY-  WORKSPACE
C     IFAIL  -INTEGER-  ERROR INDICATOR
C
C     VERSION U.K.C.  21/11/79
C     BY P.R. GRAVES-MORRIS, T.R. HOPKINS AND D.J. WINSTANLEY
C     BASED ON A PROCEDURE BY A.R. CURTIS.
C
C     EPS IS A MACHINE DEPENDENT VARIABLE.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01RAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), F(N), U(N), X(N)
      INTEGER           IW(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, W, W1, W2, W3, W4
      INTEGER           I, IFL, IP1, J, JP1, K, L, MM1, NM1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      EPS = SQRT(X02AJF())
      IFL = IFAIL
      IFAIL = 0
      IF (N.LE.0) GO TO 460
      M = 1
      NM1 = N - 1
      DO 20 I = 1, N
         IW(I) = 1
         A(I) = 0.0D0
         U(I) = X(I)
   20 CONTINUE
C     TEST FOR EQUALITY OF XI.
      IF (N.EQ.1) GO TO 80
      DO 60 J = 1, NM1
         JP1 = J + 1
         W = X(J)
         DO 40 I = JP1, N
            IF (ABS(X(I)-W).LT.EPS) GO TO 480
   40    CONTINUE
   60 CONTINUE
   80 A(1) = F(1)
      W = A(1)
      IF (N.EQ.1) RETURN
C     GI(X) IS THE REDUCED FUNCTION AT STAGE I.
C     STAGE 1 - CALCULATE G1(X) AT THE INTERPOLATION POINTS.
      DO 100 I = 2, N
         A(I) = F(I) - W
         IF (ABS(A(I)).LT.EPS) IW(I) = 0
  100 CONTINUE
C     SOLUTION IS LINEAR IF N=2.
      IF (N.EQ.2) GO TO 300
C     START BLOCK FOR THIELE COEFFICIENTS (AI, I=2,3,...,M) FOR
C     N.GT.2.
      DO 280 I = 2, NM1
         IP1 = I + 1
C        TEST THAT NOT ALL GI(XJ) ARE ZERO OR INFINITE.
         DO 120 J = I, N
            IF (IW(J).EQ.1) GO TO 140
  120    CONTINUE
         IF (M.NE.1) GO TO 320
         A(1) = F(1)
         GO TO 420
  140    W = U(I-1)
C        CALCULATE ALL POSSIBLE AI(XK) AND CHOOSE ONE NEAREST
C        UNITY IN MODULUS.
         W1 = A(J)/(U(J)-W)
         W2 = ABS(ABS(W1)-1.0D0)
         K = J
         IF (K.EQ.N) GO TO 180
         JP1 = J + 1
         DO 160 L = JP1, N
            IF (IW(L).NE.1) GO TO 160
            W3 = A(L)/(U(L)-W)
            W4 = ABS(ABS(W3)-1.0D0)
            IF (W4.GE.W2) GO TO 160
            K = L
            W1 = W3
            W2 = W4
  160    CONTINUE
  180    A(K) = A(I)
         A(I) = W1
         W = U(I)
         U(I) = U(K)
         U(K) = W
         IW(K) = IW(I)
         DO 260 J = IP1, N
C           IF GI(XJ) IS INFINITE, GO TO 9.
C           IF GI(XJ) IS ZERO, GO TO 10.
C           OTHERWISE , NORMAL ROUTE WITH IW(J)=1 IS TO 11.
            IF (IW(J)) 200, 220, 240
  200       IW(J) = 1
            A(J) = -1.0D0
            GO TO 260
  220       IW(J) = -1
            GO TO 260
  240       A(J) = A(I)*(U(J)-U(I-1))/A(J) - 1.0D0
            IF (ABS(A(J)).LT.EPS) IW(J) = 0
  260    CONTINUE
         M = I
  280 CONTINUE
      IF (IW(N).NE.1) GO TO 320
  300 M = N
      A(N) = A(N)/(U(N)-U(N-1))
C     WHEN SOLUTION IS LINEAR NO DENOMINATOR IS REQUIRED,
C     BUT A CHECK FOR ZERO SLOPE IS NECESSARY.
      IF (M.GT.2) GO TO 360
      IF (ABS(A(2)).LT.EPS) M = 1
      GO TO 420
C     IS ANY GI(XJ) INFINITE
  320 DO 340 J = I, N
         IF (IW(J).EQ.-1) GO TO 500
  340 CONTINUE
C     TEST FOR ZERO DENOMINATOR AT U(1),U(2),....,U(M).
  360 MM1 = M - 1
      DO 400 L = 1, M
         W = 0.0D0
         W1 = 1.0D0
         DO 380 J = 1, MM1
            W2 = W1 + A(J+1)*(U(L)-U(J))*W
            W = W1
            W1 = W2
  380    CONTINUE
         IF (ABS(W2).LT.EPS) GO TO 500
  400 CONTINUE
C     ZERO THE ELEMENTS OF ARRAY A THAT ARE NOT NECESSARY IN
C     DEFINING THE INTERPOLANT.
  420 IF (M.EQ.N) RETURN
      K = M + 1
      DO 440 J = K, N
         A(J) = 0.0D0
  440 CONTINUE
      RETURN
C     SET IFAIL INDICATORS.
  460 IFAIL = 1
      GO TO 540
  480 IFAIL = 2
      GO TO 540
  500 DO 520 J = 1, N
         A(J) = 0.0D0
  520 CONTINUE
      IFAIL = 3
  540 IFAIL = P01ABF(IFL,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
