      SUBROUTINE D03EDY(A,Z,U,F,LEV,LDA,NG,NGP,NGRIDX,NGRIDY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Z=CU+F on GRID(LEV) --- LU.X=Z
C
C     .. Scalar Arguments ..
      INTEGER           LDA, LEV, NG
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,7), F(NG), U(NG), Z(NG)
      INTEGER           NGP(12), NGRIDX(12), NGRIDY(12)
C     .. Local Scalars ..
      INTEGER           K, KEND, KP, NGX, NGY
C     .. Executable Statements ..
      NGX = NGRIDX(LEV)
      NGY = NGRIDY(LEV)
      KP = NGP(LEV) - NGX*NGY
      KEND = KP + NGX*NGY - NGX + 2
      Z(KP+1) = 0.0D0
      DO 20 K = KP + 2, KEND
         Z(K) = A(K,3)*A(K-1,6)*U(K+NGX-2)
   20 CONTINUE
      DO 40 K = KEND + 1, NGP(LEV)
         Z(K) = 0
   40 CONTINUE
      DO 60 K = KP + NGX, KP + NGX*NGY
         Z(K) = Z(K) + A(K,2)*A(K-NGX+1,5)*U(K-NGX+2)
   60 CONTINUE
      DO 80 K = KP + 1, NGP(LEV)
         Z(K) = Z(K) + F(K)
   80 CONTINUE
      RETURN
      END
