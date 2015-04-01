      SUBROUTINE F01BDF(N,A,IA,B,IB,DL,IFAIL)
C     MARK 3 RELEASE NAG COPYRIGHT 1972
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     MAY 1ST.  1972
C     REDUC2
C     REDUCTION OF THE GENERAL SYMMETRIC EIGENVALUE PROBLEMS
C     A*B*X=LAMBDA*X,  YT*A*B=YT*LAMBDA,
C     B*A*Y=LAMBDA*Y,  XT*B*A=XT*LAMBDA,
C     WITH SYMMETRIC MATRIX A AND SYMMETRIC POSITIVE DEFINITE
C     MATRIX
C     B, TO THE EQUIVALENT STANDARD PROBLEM Q*Z=LAMBDA*Z.
C     THE UPPER TRIANGLE, INCLUDING DIAGONAL ELEMENTS, OF A AND B
C     ARE GIVEN IN THE ARRAYS A(N,N) AND B(N,N).
C     L (B=L*LT) IS FORMED IN THE REMAINING STRICTLY LOWER TRIANGLE
C     OF THE ARRAY B WITH ITS DIAGONAL ELEMENTS IN THE ARRAY DL(N),
C     AND THE LOWER TRIANGLE OF THE SYMMETRIC MATRIX Q (Q=LT*A*L)
C     IS FORMED IN THE LOWER TRIANGLE OF THE ARRAY A, INCLUDING THE
C     DIAGONAL ELEMENTS. HENCE THE DIAGONAL ELEMENTS OF A ARE LOST.
C     THE SUBROUTINE WILL FAIL IF B, PERHAPS ON ACCOUNT OF ROUNDING
C     ERRORS, IS NOT POSITIVE DEFINITE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BDF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), DL(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, I1, ISAVE, J, J1, K, KK
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      DO 100 I = 1, N
         I1 = I - 1
         DO 80 J = I, N
            X = B(I,J)
            IF (I1.EQ.0) GO TO 40
            DO 20 KK = 1, I1
               K = I1 - KK + 1
               X = X - B(I,K)*B(J,K)
   20       CONTINUE
   40       IF (I.NE.J) GO TO 60
            IF (X.LT.0.0D0) GO TO 320
            Y = SQRT(X)
            DL(I) = Y
            GO TO 80
   60       B(J,I) = X/Y
   80    CONTINUE
  100 CONTINUE
C     L HAS BEEN FORMED IN ARRAY B
      DO 220 I = 1, N
         I1 = I + 1
         DO 200 J = 1, I
            X = A(J,I)*DL(J)
            J1 = J + 1
            IF (J1.GT.I) GO TO 140
            DO 120 K = J1, I
               X = X + A(K,I)*B(K,J)
  120       CONTINUE
  140       IF (I1.GT.N) GO TO 180
            DO 160 K = I1, N
               X = X + A(I,K)*B(K,J)
  160       CONTINUE
  180       A(I,J) = X
  200    CONTINUE
  220 CONTINUE
C     THE LOWER TRIANGLE OF A*L HAS BEEN FORMED
C     IN THE LOWER TRIANGLE OF ARRAY A
      DO 300 I = 1, N
         Y = DL(I)
         I1 = I + 1
         DO 280 J = 1, I
            X = Y*A(I,J)
            IF (I1.GT.N) GO TO 260
            DO 240 K = I1, N
               X = X + A(K,J)*B(K,I)
  240       CONTINUE
  260       A(I,J) = X
  280    CONTINUE
  300 CONTINUE
      IFAIL = 0
      RETURN
  320 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
