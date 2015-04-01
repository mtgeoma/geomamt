      SUBROUTINE G13DCP(QQ,IK,K,W,N,W2)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      INTEGER           IK, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  QQ(IK,K), W(IK,N), W2(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  XNINV
      INTEGER           I, J, T
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     CALCULATE MEAN OF EACH SERIES
C
      XNINV = 1.0D0/DBLE(N)
      DO 20 J = 1, K
         W2(J) = 0.0D0
   20 CONTINUE
      DO 60 T = 1, N
         DO 40 J = 1, K
            W2(J) = W2(J) + W(J,T)
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, K
         W2(I) = W2(I)*XNINV
   80 CONTINUE
      DO 140 T = 1, N
         DO 120 J = 1, K
            DO 100 I = J, K
               QQ(I,J) = QQ(I,J) + (W(I,T)-W2(I))*(W(J,T)-W2(J))
  100       CONTINUE
  120    CONTINUE
  140 CONTINUE
      DO 180 J = 1, K
         DO 160 I = J, K
            QQ(I,J) = QQ(I,J)*XNINV
  160    CONTINUE
  180 CONTINUE
      END
