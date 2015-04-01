      SUBROUTINE F01AYF(N,TOL,A,IA,D,E,E2)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-907 (APR 1991).
C
C     TRED3
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC
C     MATRIX, A, STORED ROW BY ROW IN THE ARRAY A(N(N+1)/2), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N-1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). E2(I) IS SET TO EQUAL
C     E(I)**2. THE ARRAY A IS USED TO STORE SUFFICIENT INFORMATION
C     FOR THE DETAILS OF THE TRANSFORMATION TO BE RECOVERABLE IN
C     THE SUBROUTINE F01AZF
C     1ST. MARCH  1972
C
C     REVISED BY VINCE FERNANDO AT MARK 14 TO INTRODUCE SCALING INTO
C     THE GENERATION OF HOUSEHOLDER MATRICES AS PROPOSED BY
C     G.W. STEWART, INTRODUCTION TO MATRIX COMPUTATIONS, CHAPTER 7.
C     TOL IS NOW A DUMMY PARAMETER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA), D(N), E(N), E2(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, HH, SCALE, SMALL
      INTEGER           I, II, IPOS, IZ, J, K, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          DDOT, IDAMAX, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSPMV, DSPR2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      DO 140 II = 1, N
         I = N - II + 1
         L = I - 1
         IZ = (I*L)/2
         H = 0.0D0
         SCALE = 0.0D0
         IF (L.EQ.0) GO TO 40
         DO 20 K = 1, L
            IPOS = IZ + K
            F = A(IPOS)
            D(K) = F
   20    CONTINUE
C        FIND THE ELEMENT OF LARGEST ABSOLUTE VALUE IN D
         K = IDAMAX(L,D,1)
         SCALE = ABS(D(K))
C        IF D IS A NULL VECTOR THEN SKIP THE TRANSFORMATION
         IF (SCALE.GE.SMALL) GO TO 60
   40    E(I) = 0.0D0
         E2(I) = 0.0D0
         H = 0.0D0
         GO TO 120
   60    CALL DSCAL(L,1.0D0/SCALE,D,1)
         H = DDOT(L,D,1,D,1)
         E2(I) = H*SCALE*SCALE
         F = D(L)
         G = SQRT(H)
         IF (F.GE.0.0D0) G = -G
         E(I) = G*SCALE
         H = H - F*G
         D(L) = F - G
         IPOS = IZ + L
         A(IPOS) = D(L)*SCALE
C        FORM A*U
         CALL DSPMV('U',L,1.0D0/H,A,D,1,0.0D0,E,1)
C        FORM P
         F = 0.0D0
         DO 80 J = 1, L
            F = F + E(J)*D(J)
   80    CONTINUE
C        FORM K
         HH = F/(H+H)
C        FORM Q
         DO 100 J = 1, L
            E(J) = E(J) - HH*D(J)
  100    CONTINUE
C        FORM REDUCED A
         CALL DSPR2('U',L,-1.0D0,D,1,E,1,A)
  120    IPOS = IZ + I
         D(I) = A(IPOS)
         A(IPOS) = H*SCALE*SCALE
  140 CONTINUE
      RETURN
      END
