      SUBROUTINE F04AQZ(N,M,RL,D,B,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     This routine originally called F04AQF.
C     LDLTSOL
C     SOLVES LDUX=B WHERE D IS A DIAGONAL MATRIX AND U IS THE
C     TRANSPOSE OF THE UNIT LOWER TRIANGULAR MATRIX L. L IS STORED
C     BY ROWS WITH THE DIAGONAL ELEMENTS OMITTED IN THE
C     ARRAY RL(I), I=1,M (M=N(N-1)/2). THE MATRIX D OCCUPIES THE N
C     ELEMENTS OF THE ARRAY D(I), I=1,N. THE SOLUTION AND RIGHT
C     HAND SIDE VECTORS ARE STORED IN X(I) AND B(I) RESPECTIVELY,
C     WHERE I=1,N.
C     1ST NOVEMBER  1972
C
C     .. Scalar Arguments ..
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), D(N), RL(M), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, II, IR, IS, IT, K, KK
C     .. Executable Statements ..
      IR = 1
      DO 60 I = 1, N
         SUM = B(I)
         IT = I - 1
         IF (IT.LT.1) GO TO 40
         DO 20 K = 1, IT
            SUM = SUM - X(K)*RL(IR)
            IR = IR + 1
   20    CONTINUE
   40    X(I) = SUM
   60 CONTINUE
      I = N + 1
      DO 120 II = 1, N
         I = I - 1
         IS = IR
         IR = IR - 1
         IT = I + 1
         SUM = X(I)/D(I)
         IF (IT.GT.N) GO TO 100
         K = N + 1
         DO 80 KK = IT, N
            K = K - 1
            SUM = SUM - X(K)*RL(IS)
            IS = IS + 2 - K
   80    CONTINUE
  100    X(I) = SUM
  120 CONTINUE
      RETURN
      END
