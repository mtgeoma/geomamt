      SUBROUTINE E02GBS(K,N,ZZ,IZR,IW,DD,V,OP,W)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     Z  IS AN  N BY N  MATRIX, AND  D  CONTAINS
C     THE NONZERO ENTRIES OF A DIAGONAL MATRIX.
C     (ZZ(2))*(D)*(ZZ(2)-TRANSP)  IS A PROJECTOR ON THE
C     SUBSPACE SPANNED BY  ZZ(2)  (THE LAST  N-K
C     COLUMNS OF  Z).  GIVEN A VECTOR  V,  THE
C     PROJECTION   (OP)=(ZZ(2))*(D)*(ZZ(2)-TRANSP)*(V)
C     IS FORMED.
C
C     W  IS A SCRATCH VECTOR.
C     ***************
C     .. Scalar Arguments ..
      INTEGER           IW, IZR, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DD(N), OP(N), V(N), W(IW), ZZ(IZR,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
      INTEGER           I, KP1, NMK
C     .. External Functions ..
      DOUBLE PRECISION  E02GBJ
      EXTERNAL          E02GBJ
C     .. Data statements ..
      DATA              ZERO/0.0D+00/
C     .. Executable Statements ..
      NMK = N - K
      KP1 = K + 1
      IF (N.LT.1) RETURN
      IF (K.GT.0) GO TO 40
C
C     ***************
C     CASE 1 ... ZZ(2)=ZZ  (K=0)
C     ***************
C
      DO 20 I = 1, N
         OP(I) = V(I)
   20 CONTINUE
      RETURN
   40 CONTINUE
      IF (K.LT.N) GO TO 80
C
C     ***************
C     CASE 2 ... ZZ(2) IS VACUOUS  (K=N)
C     ***************
C
      DO 60 I = 1, N
         OP(I) = ZERO
   60 CONTINUE
      RETURN
   80 CONTINUE
C
C     ***************
C     CASE 3 ... ZZ(2)  IS INTERMEDIATE BETWEEN
C     THE OTHER TWO CASES  (0 .LT. K .LT. N)
C     ***************
C
      DO 100 I = KP1, N
         W(I) = E02GBJ(N,ZZ(1,I),1,V,1,N,N)*DD(I)
  100 CONTINUE
      DO 120 I = 1, N
         OP(I) = E02GBJ(NMK,ZZ(I,KP1),IZR,W(KP1),1,(NMK-1)*IZR+1,NMK)
  120 CONTINUE
      RETURN
      END
