      SUBROUTINE D03EDV(A,U,V,LEV,LDA,NG,NGP,NGRIDX,NGRIDY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Solution of LU*U=V on GRID(LEV)
C     LU stored in A. V=Z from D03EDY (CTUPF) on lower levels!
C
C     .. Scalar Arguments ..
      INTEGER           LDA, LEV, NG
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,7), U(NG), V(NG)
      INTEGER           NGP(12), NGRIDX(12), NGRIDY(12)
C     .. Local Scalars ..
      INTEGER           K, KK, M
C     .. Executable Statements ..
      M = NGRIDX(LEV)
      KK = NGP(LEV) - M*NGRIDY(LEV)
      U(KK+1) = V(KK+1)
      DO 20 K = KK + 2, KK + M - 1
         U(K) = V(K) - A(K,3)*U(K-1)
   20 CONTINUE
      U(KK+M) = V(KK+M) - A(KK+M,3)*U(KK+M-1) - A(KK+M,2)*U(KK+1)
      DO 40 K = KK + M + 1, NGP(LEV)
         U(K) = V(K) - A(K,3)*U(K-1) - A(K,2)*U(K-M+1) - A(K,1)*U(K-M)
   40 CONTINUE
      KK = NGP(LEV)
      U(KK) = U(KK)/A(KK,4)
      DO 60 K = KK - 1, KK - M + 2, -1
         U(K) = (U(K)-A(K,5)*U(K+1))/A(K,4)
   60 CONTINUE
      K = KK - M + 1
      U(K) = (U(K)-A(K,5)*U(K+1)-A(K,6)*U(KK))/A(K,4)
      DO 80 K = KK - M, KK - M*NGRIDY(LEV) + 1, -1
         U(K) = (U(K)-A(K,5)*U(K+1)-A(K,6)*U(K+M-1)-A(K,7)*U(K+M))/A(K,
     *          4)
   80 CONTINUE
      RETURN
      END
