      SUBROUTINE E01AAF(A,B,C,N1,N2,N,X)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           N, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N1), B(N1), C(N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  RAI, RAK, RBI, RBK
      INTEGER           I, II, K, L, M, NN1
C     .. Executable Statements ..
      NN1 = N + 1
      DO 20 I = 1, NN1
         A(I) = A(I) - X
   20 CONTINUE
      M = 2
      K = 1
      L = 0
   40 CONTINUE
      RAK = A(K)
      RBK = B(K)
      DO 60 I = M, NN1
         II = I + L - 1
         RAI = A(I)
         RBI = (RBK*RAI-B(I)*RAK)/(RAI-RAK)
         B(I) = RBI
         C(II) = RBI
   60 CONTINUE
      L = L + NN1 - M
      K = M
      M = M + 1
      IF (NN1-K) 80, 80, 40
   80 RETURN
      END
