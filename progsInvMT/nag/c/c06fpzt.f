      SUBROUTINE C06FPZ(N,NQ,Q)
CVD$R NOVECTOR
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER           N, NQ
C     .. Array Arguments ..
      INTEGER           Q(30)
C     .. Local Scalars ..
      INTEGER           I, K, L, NN
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      NN = N
      K = 0
C
C     Trap the special case N = 1
C
      IF (N.EQ.1) THEN
         NQ = 1
         Q(1) = 1
         RETURN
      END IF
C
C     Factors of 6 --
C
   20 IF (MOD(NN,6).NE.0) GO TO 40
      K = K + 1
      Q(K) = 6
      NN = NN/6
      IF (NN.EQ.1) GO TO 160
      GO TO 20
C
C     Factors of 4 --
C
   40 IF (MOD(NN,4).NE.0) GO TO 60
      K = K + 1
      Q(K) = 4
      NN = NN/4
      IF (NN.EQ.1) GO TO 160
      GO TO 40
C
C     Factors of 2 --
C
   60 IF (MOD(NN,2).NE.0) GO TO 80
      K = K + 1
      Q(K) = 2
      NN = NN/2
      IF (NN.EQ.1) GO TO 160
      GO TO 60
C
C     Factors of 3 --
C
   80 IF (MOD(NN,3).NE.0) GO TO 100
      K = K + 1
      Q(K) = 3
      NN = NN/3
      IF (NN.EQ.1) GO TO 160
      GO TO 80
C
C     Remaining odd factors --
C
  100 L = 5
      I = 2
C
C     I is alternatively 2 or 4 --
C
  120 IF (MOD(NN,L).NE.0) GO TO 140
      K = K + 1
      Q(K) = L
      NN = NN/L
      IF (NN.EQ.1) GO TO 160
      GO TO 120
  140 L = L + I
      I = 6 - I
      GO TO 120
  160 NQ = K
C
      RETURN
      END
