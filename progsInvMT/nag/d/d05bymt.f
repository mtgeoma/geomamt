      SUBROUTINE D05BYM(V,B,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     -------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     This routine solves a Vandermonde matrix with elements
C     V(j),  j = 1,2,...,N     and  right hand B, using the
C     Bjorck & Pereyra  algorithm. The solution overwrites  B.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), V(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XK
      INTEGER           I, K
C     .. Executable Statements ..
C
      DO 40 K = 1, N - 1
         XK = V(K)
C
         DO 20 I = N, K + 1, -1
            B(I) = B(I) - XK*B(I-1)
   20    CONTINUE
   40 CONTINUE
C
      DO 100 K = N - 1, 1, -1
         DO 60 I = K + 1, N
            B(I) = B(I)/(V(I)-V(I-K))
   60    CONTINUE
C
         DO 80 I = K, N - 1
            B(I) = B(I) - B(I+1)
   80    CONTINUE
  100 CONTINUE
C
      RETURN
      END
