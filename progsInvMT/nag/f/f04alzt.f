      SUBROUTINE F04ALZ(N,M,IR,RL,IL,M1,B,IB,X,IX)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C     Originally called F04ALF
C     CHOBANDSOL
C     SOLVES AX = B, WHERE A IS A POSITIVE DEFINITE SYMMETRIC BAND
C     MATRIX WITH M LINES ON EITHER SIDE OF THE DIAGONAL AND B IS
C     AN N*IR MATRIX OF IR RIGHT-HAND SIDES. THE SUBROUTINE F04ALZ
C     MUST BE PRECEDED BY F03AGZ IN WHICH L IS PRODUCED IN L(I,K),
C     FROM A. AX=B IS SOLVED IN TWO STEPS, LY = B AND UX = Y.
C     THE MATRIX B IS RETAINED IN ORDER TO FACILITATE THE
C     REFINEMENT OF X, BUT X IS OVER WRITTEN ON Y. HOWEVER, X AND B
C     CAN BE IDENTIFIED IN THE CALL OF THE SUBROUTINE.
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IB, IL, IR, IX, M, M1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IB,IR), RL(IL,M1), X(IX,IR)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y
      INTEGER           I, II, IP, IQ, IS, J, K, KK
C     .. Executable Statements ..
      IS = M
      DO 140 J = 1, IR
C        SOLUTION OF LY = B
         DO 60 I = 1, N
            IP = 1
            IF (I.LE.(M+1)) IP = M - I + 2
            IQ = I
            Y = B(I,J)
            K = IS + 1
            IF (IS.LT.IP) GO TO 40
            DO 20 KK = IP, IS
               K = K - 1
               IQ = IQ - 1
               Y = Y - RL(I,K)*X(IQ,J)
   20       CONTINUE
   40       X(I,J) = Y*RL(I,M+1)
   60    CONTINUE
C        SOLUTION OF UX = Y
         I = N + 1
         DO 120 II = 1, N
            I = I - 1
            IP = 1
            IF ((N-I).LE.M) IP = M - N + I + 1
            Y = X(I,J)
            IQ = I
            K = IS + 1
            IF (IS.LT.IP) GO TO 100
            DO 80 KK = IP, IS
               K = K - 1
               IQ = IQ + 1
               Y = Y - RL(IQ,K)*X(IQ,J)
   80       CONTINUE
  100       X(I,J) = Y*RL(I,M+1)
  120    CONTINUE
  140 CONTINUE
      RETURN
      END
