      SUBROUTINE F01AGF(N,ATOL,A,IA,D,E,E2)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-905 (APR 1991).
C
C     TRED1
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N - 1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). E2(I) IS SET TO EQUAL
C     E(I)**2. THE LOWER TRIANGLE OF THE ARRAY A, TOGETHER WITH THE
C     ARRAY E, IS USED TO STORE SUFFICIENT INFORMATION FOR THE
C     DETAILS OF THE TRANSFORMATION TO BE RECOVERABLE IN THE
C     SUBROUTINE F01AHF. THE STRICTLY UPPER TRIANGLE OF THE ARRAY A
C     IS LEFT UNALTERED.
C     1ST AUGUST 1971
C
C     REVISED BY VINCE FERNANDO AT MARK 14 TO INTRODUCE SCALING INTO
C     THE GENERATION OF HOUSEHOLDER MATRICES AS PROPOSED BY
C     G.W. STEWART, INTRODUCTION TO MATRIX COMPUTATIONS, CHAPTER 7.
C     ATOL IS NOW A DUMMY PARAMETER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ATOL
      INTEGER           IA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N), E(N), E2(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, SCALE, SMALL
      INTEGER           I, II, J, K, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          DDOT, IDAMAX, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSYMV, DSYR2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      DO 20 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
   20 CONTINUE
      DO 160 II = 1, N
         I = N - II + 1
         L = I - 1
         H = 0.0D0
         IF (L.EQ.0) GO TO 60
C        FIND THE ELEMENT OF LARGEST ABSOLUTE VALUE IN D
         K = IDAMAX(L,D,1)
         SCALE = ABS(D(K))
C        IF D IS A NULL VECTOR THEN SKIP THE TRANSFORMATION
         IF (SCALE.GE.SMALL) GO TO 80
         DO 40 J = 1, L
            D(J) = A(I-1,J)
            A(I-1,J) = A(I,J)
            A(I,J) = 0.0D0
   40    CONTINUE
   60    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 160
   80    CALL DSCAL(L,1.0D0/SCALE,D,1)
         H = DDOT(L,D,1,D,1)
         E2(I) = H*SCALE*SCALE
         F = D(I-1)
         G = SQRT(H)
         IF (F.GE.0.0D0) G = -G
         E(I) = G*SCALE
         H = H - F*G
         D(I-1) = F - G
C        FORM A*U
         CALL DSYMV('L',L,1.0D0/H,A,IA,D,1,0.0D0,E,1)
C        FORM P
         F = 0.0D0
         DO 100 J = 1, L
            F = F + E(J)*D(J)
  100    CONTINUE
C        FORM K
         H = F/(H+H)
C        FORM Q
         DO 120 J = 1, L
            E(J) = E(J) - H*D(J)
  120    CONTINUE
C        FORM REDUCED A
         CALL DSYR2('L',L,-1.0D0,D,1,E,1,A,IA)
         DO 140 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F*SCALE
  140    CONTINUE
  160 CONTINUE
      RETURN
      END
