      SUBROUTINE F02AMF(N,EPS,D,E,Z,IZ,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     TQL2
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS OF A
C     TRIDIAGONAL MATRIX, T, GIVEN WITH ITS DIAGONAL ELEMENTS IN
C     THE ARRAY D(N) AND ITS SUB-DIAGONAL ELEMENTS IN THE LAST N
C     - 1 STORES OF THE ARRAY E(N), USING QL TRANSFORMATIONS. THE
C     EIGENVALUES ARE OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE
C     ARRAY D IN ASCENDING ORDER. THE EIGENVECTORS ARE FORMED IN
C     THE ARRAY Z(N,N), OVERWRITING THE ACCUMULATED
C     TRANSFORMATIONS AS SUPPLIED BY THE SUBROUTINE F01AJF. THE
C     SUBROUTINE WILL FAIL IF ALL EIGENVALUES TAKE MORE THAN 30*N
C     ITERATIONS.
C     1ST APRIL 1972
C
C     .. Parameters ..
      INTEGER           VLEN
      PARAMETER         (VLEN=128)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AMF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IFAIL, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), E(N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, F, G, H, P, R, S
      INTEGER           I, I1, IPOS, ISAVE, ISEG, J, K, L, M
C     .. Local Arrays ..
      DOUBLE PRECISION  CC(VLEN), SS(VLEN)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06QXF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
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
      DO 300 L = 1, N
         H = EPS*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C        LOOK FOR SMALL SUB-DIAG ELEMENT
         DO 60 M = L, N
            IF (ABS(E(M)).LE.B) GO TO 80
   60    CONTINUE
   80    IF (M.EQ.L) GO TO 280
  100    IF (J.LE.0) GO TO 400
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
         DO 260 K = M - 1, L, -VLEN
            ISEG = MAX(K-VLEN+1,L)
            DO 240 I = K, ISEG, -1
               G = C*E(I)
               H = C*P
               IF (ABS(P).LT.ABS(E(I))) GO TO 200
               C = E(I)/P
               R = SQRT(C*C+1.0D0)
               E(I+1) = S*P*R
               S = C/R
               C = 1.0D0/R
               GO TO 220
  200          C = P/E(I)
               R = SQRT(C*C+1.0D0)
               E(I+1) = S*E(I)*R
               S = 1.0D0/R
               C = C/R
  220          P = C*D(I) - S*G
               D(I+1) = H + S*(C*G+S*D(I))
C           STORE ROTATIONS
               CC(VLEN-K+I) = C
               SS(VLEN-K+I) = -S
  240       CONTINUE
C        UPDATE VECTORS
            IPOS = VLEN - K + ISEG
            CALL F06QXF('Right','Variable','Backward',N,K-ISEG+2,1,
     *                  K-ISEG+2,CC(IPOS),SS(IPOS),Z(1,ISEG),IZ)
  260    CONTINUE
         E(L) = S*P
         D(L) = C*P
         IF (ABS(E(L)).GT.B) GO TO 100
  280    D(L) = D(L) + F
  300 CONTINUE
C     ORDER EIGENVALUES AND EIGENVECTORS
      DO 380 I = 1, N
         K = I
         P = D(I)
         I1 = I + 1
         IF (I1.GT.N) GO TO 340
         DO 320 J = I1, N
            IF (D(J).GE.P) GO TO 320
            K = J
            P = D(J)
  320    CONTINUE
  340    IF (K.EQ.I) GO TO 380
         D(K) = D(I)
         D(I) = P
         DO 360 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  360    CONTINUE
  380 CONTINUE
      IFAIL = 0
      RETURN
  400 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
