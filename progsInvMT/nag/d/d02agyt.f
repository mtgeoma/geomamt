      SUBROUTINE D02AGY(Y,X,H,N,N1,AUX,WSPACE,WSPAC1,PARAM)
C
C     USES THE TECHNIQUE OF NAG LIBRARY
C     PROCEDURE D02AAF WITH THE
C     SPECIFICATION OF AUX CHANGED
C     ALL IMPLICITLY DECLARED REALS MAY
C     BE DECLARED DOUBLE PRECISION
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, X
      INTEGER           N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  PARAM(N1), WSPAC1(N), WSPACE(N,9), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          AUX
C     .. Local Scalars ..
      DOUBLE PRECISION  C1, C2, DUM, U, V, W
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      CALL AUX(WSPAC1,Y,X,PARAM)
      C1 = 1.0D0/3.0D0
      C2 = 1.0D0/6.0D0
      DO 20 I = 1, N
         WSPACE(I,3) = WSPAC1(I)
   20 CONTINUE
      U = C1*H
      DO 40 I = 1, N
         WSPACE(I,4) = Y(I)
         Y(I) = Y(I) + U*WSPACE(I,3)
   40 CONTINUE
      CALL AUX(WSPAC1,Y,U+X,PARAM)
      DO 60 I = 1, N
         WSPACE(I,5) = WSPAC1(I)
   60 CONTINUE
      V = H*C2
      DO 80 I = 1, N
         Y(I) = WSPACE(I,4) + V*(WSPACE(I,3)+WSPACE(I,5))
   80 CONTINUE
      CALL AUX(WSPAC1,Y,U+X,PARAM)
      DO 100 I = 1, N
         WSPACE(I,5) = WSPAC1(I)
  100 CONTINUE
      U = H*0.125D0
      V = H*0.375D0
      DO 120 I = 1, N
         Y(I) = WSPACE(I,4) + WSPACE(I,3)*U + WSPACE(I,5)*V
  120 CONTINUE
      U = 0.5D0*H
      V = 1.5D0*H
      W = 2.0D0*H
      CALL AUX(WSPAC1,Y,U+X,PARAM)
      DO 140 I = 1, N
         Y(I) = WSPACE(I,4) + WSPACE(I,3)*U - WSPACE(I,5)*V + WSPAC1(I)
     *          *W
         WSPACE(I,5) = WSPAC1(I)
  140 CONTINUE
      X = X + H
      CALL AUX(WSPAC1,Y,X,PARAM)
      U = H*C2
      V = 2.0D0*H*C1
      DO 160 I = 1, N
         W = WSPACE(I,4) + U*(WSPACE(I,3)+WSPAC1(I)) + WSPACE(I,5)*V
         DUM = W - Y(I)
         WSPACE(I,2) = 0.2D0*ABS(DUM)
         Y(I) = W
  160 CONTINUE
      RETURN
      END
