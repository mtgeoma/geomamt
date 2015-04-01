      SUBROUTINE G13DSV(K,I,ST,L1,L3,ITW,M,ACFVAR,IM,IK,C,LP,MD)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           I, IK, IM, ITW, K, L1, L3, LP, M, MD
C     .. Array Arguments ..
      DOUBLE PRECISION  ACFVAR(IM,MD*K*K), C(IK,IK,MD)
      CHARACTER*1       ST(80)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, SUM2
      INTEGER           J, L, L2, LLP
      CHARACTER*81      REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Executable Statements ..
C
C     This subroutine will print out the significant cross-correlations
C     between lags L1 and L3
C
      DO 20 J = 1, 80
         ST(J) = ' '
   20 CONTINUE
      DO 60 J = 1, K
         ST((J-1)*(M+3)+1) = '*'
         LLP = 0
         DO 40 L = L1, L3
            LLP = LLP + 1
            L2 = (L-1)*K*K + (J-1)*K + I
            SUM = 1.96D0*ACFVAR(L2,L2)
            L2 = (J-1)*(M+3) + LLP + 2
            SUM2 = C(I,J,L)
            IF (SUM2.GT.SUM) THEN
               ST(L2) = '+'
            ELSE IF (SUM2.LT.-SUM) THEN
               ST(L2) = '-'
            ELSE
               ST(L2) = '.'
            END IF
   40    CONTINUE
   60 CONTINUE
      ST(ITW) = '*'
      WRITE (REC,FMT=99999) (ST(J),J=1,80)
      CALL X04BAF(LP,REC)
      RETURN
C
99999 FORMAT (' ',80A1)
      END
