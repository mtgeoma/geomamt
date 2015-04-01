      SUBROUTINE F04MAY(N,A,INJ,IAJ,D,IK,B,LROW)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS SUBROUTINE COMPUTES THE PRODUCT VECTOR OF THE INVERSE OF THE
C     PRECONDITIONING MATRIX AND THE VECTOR B. THE RESULT IS PLACED
C     IN B.
C
C     .. Scalar Arguments ..
      INTEGER           IAJ, LROW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IAJ), B(N), D(N)
      INTEGER           IK(N,2), INJ(IAJ)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIC, BIR
      INTEGER           IC, IIP, IPI, IR, K, KL, KP
C     .. Executable Statements ..
      KP = 1
C
C     PERFORM THE FORWARD SUBSTITUTION.
C
      DO 60 IIP = 1, N
         IC = IK(IIP,2)
         KL = KP + IK(IC,1) - 1
         BIC = B(IC)
         IF (KP.GT.KL) GO TO 40
         DO 20 K = KP, KL
            IR = INJ(K)
            B(IR) = B(IR) - A(K)*BIC
   20    CONTINUE
   40    KP = KL + 1
   60 CONTINUE
      KL = LROW
C
C     PERFORM THE BACK SUBSTITUTION.
C
      DO 120 IPI = 1, N
         IIP = N + 1 - IPI
         IR = IK(IIP,2)
         BIR = 0.0D0
         KP = KL - IK(IR,1) + 1
         IF (KP.GT.KL) GO TO 100
         DO 80 K = KP, KL
            IC = INJ(K)
            BIR = BIR - A(K)*B(IC)
   80    CONTINUE
  100    B(IR) = B(IR)/D(IR) + BIR
         KL = KP - 1
  120 CONTINUE
      RETURN
C
C     END OF F04MAY. (MA31G )
C
      END
