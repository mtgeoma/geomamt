      SUBROUTINE G13DCS(IFLAG,K7,X,F,G,NSTATE,IW,W2)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     DUMMY ROUTINE OBJFUN REQUIRED BY E04XAF WHICH CALLS FUNCT
C     OF E04JBL
C
C     .. Scalars in Common ..
C
      INTEGER           LIW, LW, N4
C     .. Scalar Arguments ..
      DOUBLE PRECISION  F
      INTEGER           IFLAG, K7, NSTATE
C     .. Array Arguments ..
      DOUBLE PRECISION  G(K7), W2(LW), X(K7)
      INTEGER           IW(LIW)
C     .. Local Scalars ..
      INTEGER           I, I2, LEW10, LEW11, LX1, LX2, N2
C     .. External Subroutines ..
      EXTERNAL          G13DCR
C     .. Common blocks ..
      COMMON            /DG13DC/N4, LIW, LW
C     .. Executable Statements ..
      N2 = IW(1)
      LEW10 = IW(2)
      LEW11 = LEW10 + N2
      LX1 = 2*N2 + 1
      LX2 = LX1 + N2
C
C     SET 'X' ARRAY FOR G13DCR INTO W2(LX1) (AND GRADIENT VECTOR
C     INTO W2(LX2)
C
      DO 20 I = 1, N2
         W2(LX1-1+I) = W2(LEW10-1+I)
   20 CONTINUE
C
      I2 = 1
      DO 40 I = 1, N4
         IF (W2(4*N2+I).LT.1.0D0) GO TO 40
         W2(LX1-1+I) = X(I2)
         I2 = I2 + 1
   40 CONTINUE
C
      CALL G13DCR(IFLAG,N2,W2(LX1),F,W2(LX2),IW,LIW,W2,LW)
C
      RETURN
      END
