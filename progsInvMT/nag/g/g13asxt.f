      SUBROUTINE G13ASX(PAR,P,Q,PS,QS,FJAC,VAR,ACFVAR,NPAR,M,IM,ISEA,W,
     *                  N,LMAX,IFAILA)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IFAILA, IM, ISEA, LMAX, M, N, NPAR, P, PS, Q, QS
C     .. Array Arguments ..
      DOUBLE PRECISION  ACFVAR(IM,M), FJAC(M,NPAR), PAR(NPAR),
     *                  VAR(NPAR+1,NPAR), W(LMAX*ISEA+M)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, SUM2
      INTEGER           I, J, J2, K, LW1
C     .. External Subroutines ..
      EXTERNAL          F01ADF, G13ASY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      LW1 = LMAX*ISEA + 1
C
C     set up the 4 sets of columns of FJAC
C
      IF (P.GT.0) CALL G13ASY(PAR,P,FJAC,1,0,W,M,NPAR,W(LW1),P)
      IF (Q.GT.0) CALL G13ASY(PAR,Q,FJAC,1,P,W,M,NPAR,W(LW1),Q)
      IF (PS.GT.0) CALL G13ASY(PAR,PS,FJAC,ISEA,P+Q,W,M,NPAR,W(LW1),
     *                         PS*ISEA)
      IF (QS.GT.0) CALL G13ASY(PAR,QS,FJAC,ISEA,P+Q+PS,W,M,NPAR,W(LW1),
     *                         QS*ISEA)
C
      DO 40 J = 1, P
         DO 20 I = 1, M
            FJAC(I,J) = -FJAC(I,J)
   20    CONTINUE
   40 CONTINUE
C
      DO 80 J = 1, PS
         DO 60 I = 1, M
            FJAC(I,J+P+Q) = -FJAC(I,J+P+Q)
   60    CONTINUE
   80 CONTINUE
C
C     set upper triangle of VAR to inverse of (FJAC' * FJAC)
C
      DO 140 I = 1, NPAR
         DO 120 J = 1, I
            SUM = 0.0D0
            DO 100 K = 1, M
               SUM = SUM + FJAC(K,I)*FJAC(K,J)
  100       CONTINUE
            VAR(J,I) = SUM
  120    CONTINUE
  140 CONTINUE
C
      IFAILA = 1
      CALL F01ADF(NPAR,VAR,NPAR+1,IFAILA)
      IF (IFAILA.GT.0) RETURN
C
C     set VAR to inverse
C
      DO 180 I = 1, NPAR
         DO 160 J = 1, I
            VAR(J,I) = VAR(I+1,J)
  160    CONTINUE
  180 CONTINUE
C
      DO 220 I = 1, NPAR
         DO 200 J = 1, I
            VAR(I,J) = VAR(J,I)
  200    CONTINUE
  220 CONTINUE
C
      SUM2 = 1.0D0/DBLE(N)
      DO 300 I = 1, M
         DO 280 J2 = 1, I
            SUM = 0.0D0
            DO 260 K = 1, NPAR
               DO 240 J = 1, NPAR
                  SUM = SUM + FJAC(I,K)*VAR(K,J)*FJAC(J2,J)
  240          CONTINUE
  260       CONTINUE
            IF (I.EQ.J2) THEN
               SUM = 1.0D0 - SUM
            ELSE
               SUM = -SUM
            END IF
            ACFVAR(I,J2) = SUM*SUM2
            ACFVAR(J2,I) = SUM*SUM2
C
  280    CONTINUE
  300 CONTINUE
      RETURN
C
      END
