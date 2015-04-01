      SUBROUTINE D05BYH(LENP,IORDER,LIQ,WT,P,R)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine finds the square root weights
C     of BDF reducible quadrature of orders 4 to 6.
C     On exit, WT contains 2**(LENP+1) square
C     root weights of the method.
C     ++++++++++++++++++++++++++++++++++++++++++++++++
C     ------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IORDER, LENP, LIQ
C     .. Array Arguments ..
      DOUBLE PRECISION  P(0:2*LIQ-1), R(0:2*LIQ-1), WT(0:2*LIQ-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, J, M, NLEN, NLENM1, NN, NNM1, NNN, NNNM1
C     .. Local Arrays ..
      DOUBLE PRECISION  ALFA(0:6)
C     .. External Subroutines ..
      EXTERNAL          D05BYN, D05BYP
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, SQRT
C     .. Executable Statements ..
C
      CALL D05BYN(IORDER,ALFA)
C
      WT(0) = 1.D0/SQRT(ALFA(0))
      WT(1) = -0.5D0*ALFA(1)*(WT(0)**3)
C
      DO 20 I = 0, IORDER
         ALFA(I) = -0.5D0*ALFA(I)
   20 CONTINUE
C
      DO 120 M = 1, LENP
         NLEN = 2**M
         NLENM1 = NLEN - 1
         NN = 2*NLEN
         NNM1 = NN - 1
         NNN = 2*NN
         NNNM1 = NNN - 1
C
         DO 40 I = NLEN, NNNM1
            WT(I) = 0.D0
   40    CONTINUE
C
         DO 60 I = 0, NNNM1
            P(I) = WT(I)
   60    CONTINUE
C
         CALL D05BYP(R,P,NNN)
C
         DO 100 I = NLEN, NN
            SUM = 0.0D0
            NNN = MIN(IORDER,I)
            DO 80 J = 0, NNN
               SUM = SUM + R(I-J)*ALFA(J)
   80       CONTINUE
            WT(I) = SUM
  100    CONTINUE
C
  120 CONTINUE
C
      RETURN
      END
