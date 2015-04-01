      SUBROUTINE G02AAS(UPLO,N,A,B,C)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes C = AB where A, B and C are symmetric matrices in
C     packed storage (lower triangular by columns)
C
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*), C(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, II, IJ, J, JJ, K, KI, KJ
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. Executable Statements ..
      IF (N.EQ.1) THEN
         C(1) = A(1)*B(1)
      ELSE IF (N.GT.1) THEN
         IF (UPLO.EQ.'L' .OR. UPLO.EQ.'l') THEN
            IJ = 1
            DO 80 I = 1, N
               JJ = 0
               DO 60 J = I, N
                  SUM = 0.0D0
                  KI = I
                  KJ = J
                  DO 20 K = 1, I - 1
                     WRITE (6,FMT=*) KI, KJ
                     SUM = SUM + A(KJ)*B(KI)
                     KI = KI + N - K
                     KJ = KJ + N - K
   20             CONTINUE
                  DO 40 K = I, J - 1
                     WRITE (6,FMT=*) KI, KJ
                     SUM = SUM + A(KJ)*B(KI)
                     KI = KI + 1
                     KJ = KJ + N - K
   40             CONTINUE
                  WRITE (6,FMT=*) KI, KJ, ' *'
                  SUM = SUM + DDOT(N-J+1,A(KJ),1,B(KI),1)
                  C(IJ) = SUM
                  IJ = IJ + 1
   60          CONTINUE
   80       CONTINUE
         ELSE IF (UPLO.EQ.'U' .OR. UPLO.EQ.'u') THEN
            IJ = 1
            II = 0
            DO 160 I = 1, N
               JJ = 0
               DO 140 J = 1, I
                  SUM = DDOT(J,A(II+1),1,B(JJ+1),1)
                  KI = II + J
                  KJ = JJ + J
                  DO 100 K = J, I - 1
                     KI = KI + 1
                     KJ = KJ + K
                     SUM = SUM + A(KI)*B(KJ)
  100             CONTINUE
                  DO 120 K = I, N - 1
                     KI = KI + K
                     KJ = KJ + K
                     SUM = SUM + A(KI)*B(KJ)
  120             CONTINUE
                  C(IJ) = SUM
                  IJ = IJ + 1
                  JJ = JJ + J
  140          CONTINUE
               II = II + I
  160       CONTINUE
         END IF
      END IF
      RETURN
      END
