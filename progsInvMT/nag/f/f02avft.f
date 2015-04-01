      SUBROUTINE F02AVF(N,ACHEPS,D,E,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 9 REVISED. IER-326 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     TQL1
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A TRIDIAGONAL
C     MATRIX,
C     T, GIVEN WITH ITS DIAGONAL ELEMENTS IN THE ARRAY D(N) AND
C     ITS SUBDIAGONAL ELEMENTS IN THE LAST N - 1 STORES OF THE
C     ARRAY E(N), USING QL TRANSFORMATIONS. THE EIGENVALUES ARE
C     OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE ARRAY D IN
C     ASCENDING ORDER. THE SUBROUTINE WILL FAIL IF ALL
C     EIGENVALUES TAKE MORE THAN 30*N ITERATIONS.
C     1ST APRIL 1972
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AVF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACHEPS
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), E(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, F, G, H, P, R, S
      INTEGER           I, I1, II, ISAVE, J, L, M, M1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I = 2, N
         E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      J = 30*N
      DO 340 L = 1, N
         H = ACHEPS*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C        LOOK FOR SMALL SUB DIAGONAL ELEMENT
         DO 60 M = L, N
            IF (ABS(E(M)).LE.B) GO TO 80
   60    CONTINUE
   80    IF (M.EQ.L) GO TO 260
  100    IF (J.LE.0) GO TO 360
         J = J - 1
C        FORM SHIFT
         G = D(L)
         H = D(L+1) - G
         IF (ABS(H).GE.ABS(E(L))) GO TO 120
         P = H*0.5D0/E(L)
         R = SQRT(P*P+1.0D0)
         H = P + R
         IF (P.LT.0.0D0) H = P - R
         D(L) = E(L)/H
         GO TO 140
  120    P = 2.0D0*E(L)/H
         R = SQRT(P*P+1.0D0)
         D(L) = E(L)*P/(1.0D0+R)
  140    H = G - D(L)
         I1 = L + 1
         IF (I1.GT.N) GO TO 180
         DO 160 I = I1, N
            D(I) = D(I) - H
  160    CONTINUE
  180    F = F + H
C        QL TRANSFORMATION
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         M1 = M - 1
         DO 240 II = L, M1
            I = M1 - II + L
            G = C*E(I)
            H = C*P
            IF (ABS(P).LT.ABS(E(I))) GO TO 200
            C = E(I)/P
            R = SQRT(C*C+1.0D0)
            E(I+1) = S*P*R
            S = C/R
            C = 1.0D0/R
            GO TO 220
  200       C = P/E(I)
            R = SQRT(C*C+1.0D0)
            E(I+1) = S*E(I)*R
            S = 1.0D0/R
            C = C/R
  220       P = C*D(I) - S*G
            D(I+1) = H + S*(C*G+S*D(I))
  240    CONTINUE
         E(L) = S*P
         D(L) = C*P
         IF (ABS(E(L)).GT.B) GO TO 100
  260    P = D(L) + F
C        ORDER EIGENVALUE
         IF (L.EQ.1) GO TO 300
         DO 280 II = 2, L
            I = L - II + 2
            IF (P.GE.D(I-1)) GO TO 320
            D(I) = D(I-1)
  280    CONTINUE
  300    I = 1
  320    D(I) = P
  340 CONTINUE
      IFAIL = 0
      RETURN
  360 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
