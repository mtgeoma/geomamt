      SUBROUTINE F01AJF(N,ATOL,A,IA,D,E,Z,IZ)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-906 (APR 1991).
C
C     TRED2
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N - 1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). THE TRANSFORMATION
C     MATRICES ARE ACCUMULATED IN THE ARRAY Z(N,N). THE ARRAY
C     A IS LEFT UNALTERED UNLESS THE ACTUAL PARAMETERS
C     CORRESPONDING TO A AND Z ARE IDENTICAL.
C     1ST AUGUST 1971
C
C     REVISED BY VINCE FERNANDO AT MARK 14 TO INTRODUCE SCALING INTO
C     THE GENERATION OF HOUSEHOLDER MATRICES AS PROPOSED BY
C     G.W. STEWART, INTRODUCTION TO MATRIX COMPUTATIONS, CHAPTER 7.
C     ATOL IS NOW A DUMMY PARAMETER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ATOL
      INTEGER           IA, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N), E(N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, HH, SCALE, SMALL
      INTEGER           I, II, J, K, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          DDOT, IDAMAX, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, DSCAL, DSYMV, DSYR2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      DO 40 I = 1, N
         DO 20 J = I, N
            Z(J,I) = A(J,I)
   20    CONTINUE
         D(I) = A(N,I)
   40 CONTINUE
      IF (N.EQ.1) GO TO 340
      DO 240 II = 2, N
         I = N - II + 2
         L = I - 1
         IF (L.EQ.1) GO TO 60
C        FIND THE ELEMENT OF LARGEST ABSOLUTE VALUE IN D
         K = IDAMAX(L,D,1)
         SCALE = ABS(D(K))
C        IF D IS A NULL VECTOR THEN SKIP THE TRANSFORMATION
         IF (SCALE.GE.SMALL) GO TO 120
   60    E(I) = D(L)
         H = 0.0D0
         DO 80 J = 1, L
            Z(J,I) = 0.0D0
   80    CONTINUE
         DO 100 J = 1, L
            Z(I,J) = 0.0D0
            D(J) = Z(I-1,J)
  100    CONTINUE
         GO TO 220
  120    CALL DSCAL(L,1.0D0/SCALE,D,1)
         H = DDOT(L,D,1,D,1)
         F = D(I-1)
         G = SQRT(H)
         IF (F.GE.0.0D0) G = -G
         E(I) = G*SCALE
         H = H - F*G
         D(I-1) = F - G
C        COPY U
         DO 140 J = 1, L
            Z(J,I) = D(J)
  140    CONTINUE
C        FORM A*U
         CALL DSYMV('L',L,1.0D0/H,Z,IZ,D,1,0.0D0,E,1)
C        FORM P
         F = 0.0D0
         DO 160 J = 1, L
            F = F + E(J)*D(J)
  160    CONTINUE
C        FORM K
         HH = F/(H+H)
C        FORM Q
         DO 180 J = 1, L
            E(J) = E(J) - HH*D(J)
  180    CONTINUE
C        FORM REDUCED A
         CALL DSYR2('L',L,-1.0D0,D,1,E,1,Z,IZ)
         DO 200 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  200    CONTINUE
  220    D(I) = H
  240 CONTINUE
C     ACCUMULATION OF TRANSFORMATION MATRICES
      DO 300 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H.EQ.0.0D0) GO TO 260
         CALL DGEMV('T',L,L,1.0D0/H,Z,IZ,Z(1,I),1,0.0D0,D,1)
         CALL DGER(L,L,-1.0D0,Z(1,I),1,D,1,Z,IZ)
  260    DO 280 J = 1, L
            Z(J,I) = 0.0D0
  280    CONTINUE
  300 CONTINUE
      DO 320 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
  320 CONTINUE
  340 Z(N,N) = 1.0D0
      E(1) = 0.0D0
      RETURN
      END
