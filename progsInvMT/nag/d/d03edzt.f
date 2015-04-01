      SUBROUTINE D03EDZ(A,Z,U,V,NGX,NGXY,LDA)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     V=C(U-V) on finest grid - simple way to find residual.
C     Z is used as scratch space
C
C     .. Scalar Arguments ..
      INTEGER           LDA, NGX, NGXY
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,7), U(NGXY), V(NGXY), Z(NGXY)
C     .. Local Scalars ..
      INTEGER           K
C     .. Executable Statements ..
      DO 20 K = 2, NGXY
         Z(K) = U(K) - V(K)
   20 CONTINUE
      V(1) = 0
      DO 40 K = 2, NGXY - NGX + 2
         V(K) = A(K,3)*A(K-1,6)*Z(K+NGX-2)
   40 CONTINUE
      DO 60 K = NGXY - NGX + 3, NGXY
         V(K) = 0
   60 CONTINUE
      DO 80 K = NGX, NGXY
         V(K) = V(K) + A(K,2)*A(K-NGX+1,5)*Z(K-NGX+2)
   80 CONTINUE
      RETURN
      END
