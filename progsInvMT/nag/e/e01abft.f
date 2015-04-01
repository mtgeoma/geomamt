      SUBROUTINE E01ABF(N,P,A,G,N1,N2,IFAIL)
C     MARK 2 RELEASE.  NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 7A REVISED  IER-158 (MAR 1979)
C     MARK 9 REVISED. IER-351 (SEP 1981)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01ABF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IFAIL, N, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N1), G(N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, D, F, PROD, RK, W, X, Y
      INTEGER           I, I1, ISAVE, J, K, L, M, NN1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      NN1 = N + 1
      ISAVE = IFAIL
      IFAIL = 0
      IF (ABS(P)-1.0D0) 40, 20, 20
   20 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      GO TO 260
   40 I = 1
      K = N + N - 1
      W = A(N)
      B = (1.0D0-P)*W
      G(1) = W
      W = A(NN1)
      B = B + W*P
      G(2) = W
   60 M = 1
   80 CONTINUE
      IF (K.LT.1) GO TO 240
      DO 100 J = 1, K
         A(J) = A(J) - A(J+1)
  100 CONTINUE
      IF (M-2) 120, 140, 120
  120 M = 2
      K = K - 1
      GO TO 80
  140 L = N - I
      C = A(L)
      D = A(L+1)
      L = I + I + 1
      G(L) = C
      G(L+1) = D
      Y = I
      RK = P + Y
      PROD = RK
      M = 1
      L = I + I
  160 I1 = L + 1
      PROD = PROD/DBLE(I1)
      DO 180 J = 1, L
         I1 = I1 - 1
         X = J
         PROD = ((RK-X)/DBLE(I1))*PROD
  180 CONTINUE
      F = PROD*D
      B = B + F
      IF (M-2) 200, 220, 200
  200 M = 2
      RK = 1.0D0 - P + Y
      PROD = RK
      D = C
      GO TO 160
  220 IF ((N-I).EQ.1) GO TO 240
      I = I + 1
      GO TO 60
  240 L = N + N + 1
      G(L) = B
  260 RETURN
      END
