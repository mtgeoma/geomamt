      SUBROUTINE G13DCR(IFLAG,N2,X,SUM,GC,IW,LIW,W2,LW)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     SUBROUTINE FUNCT OF E04JBL WHICH CALLS NFUNCT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SUM
      INTEGER           IFLAG, LIW, LW, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  GC(N2), W2(LW), X(N2)
      INTEGER           IW(LIW)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ADDLOG, COND, LMAX, NORM
      INTEGER           K, K4, K5, K6, KR, LEW6, LEW7, LIW1, LP, LW1,
     *                  LW10, LW11, LW12, LW13, LW14, LW15, LW16, LW17,
     *                  LW18, LW19, LW2, LW20, LW21, LW3, LW4, LW5, LW6,
     *                  LW7, LW8, LW9, NITER, P, Q
      LOGICAL           FULLP, FULLQ, MEAN, NOPRIN
C     .. External Subroutines ..
      EXTERNAL          G13DCY
C     .. Common blocks ..
      COMMON            /AG13DC/LW6, LW7, LW9, LW10, LW15, LW19, LW20,
     *                  LW21, LIW1
      COMMON            /BG13DC/ADDLOG, LMAX, COND, NORM, K, P, Q, K4,
     *                  K5, LW14, LW18, LP, K6, NITER, LEW6, LEW7, LW2,
     *                  LW1, LW3, LW4, LW5, LW8, LW11, LW12, LW13, LW16,
     *                  LW17, KR, MEAN, NOPRIN, FULLP, FULLQ
C     .. Executable Statements ..
      CALL G13DCY(X,N2,W2(LW1),W2(LW2),W2(LW3),W2(LW4),W2(LW5),W2(LW6),
     *            W2(LW7),W2(LW8),W2(LW9),W2(LW10),W2(LW11),W2(LW12),
     *            W2(LW13),W2(LW14),W2(LW15),W2(LW16),W2(LW17),W2(LW18),
     *            W2(LW19),W2(LW20),W2(LW21),W2(LEW6),W2(LEW7),IW(LIW1),
     *            SUM,IFLAG)
C
      RETURN
      END
