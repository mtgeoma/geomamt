      SUBROUTINE G13ASF(N,V,MR,M,PAR,NPAR,ISHOW,C,ACFVAR,IM,SUM2,IDF,
     *                  SIGLEV,INTGR,LMAX,WORK,LWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15 REVISED. IER-942 (APR 1991).
C     MARK 16 REVISED. IER-1121 (JUL 1993).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13ASF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGLEV, SUM2
      INTEGER           IDF, IFAIL, IM, ISHOW, LMAX, LWORK, M, N, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  ACFVAR(IM,M), C(M), PAR(NPAR), V(N), WORK(LWORK)
      INTEGER           INTGR(LMAX), MR(7)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGEST, LIMIT, QSTAT, SUM, XM, XV
      INTEGER           I, I2, IERR, IFAILA, IFAULT, IS, ISEA, J, K, L,
     *                  L2, L3, LP, LW2, LW3, P, PS, Q, QS
      LOGICAL           SEASON, SHOW, STAT
C     .. Local Arrays ..
      CHARACTER*80      P01REC(4), REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05HDZ, G13ABF, G13ASX, X04ABF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT
C     .. Executable Statements ..
C
C     ******************************************************************
C
C     This routine calculates the variance-covariance matrix of the
C     residual autocorrelation function in a univariate ARMA model.
C     It also returns a portmanteau statistic along with its'
C     significance level. A user would normally summon this routine
C     after calling either G13AEF, G13AFF, G13BEF or G13DCF to fit a
C     univariate ARMA model.
C
C     The computational procedure is described in a paper by
C     A.I. McLeod
C       'On the distribution of residual autocorrelations in Box-Jenkins
C        models'
C        Journal of the Royal Statistical Society, Series B, Vol 40,
C        pages 296 - 302, 1978.
C
C     ******************************************************************
C
C     first test for errors in the input arguments
C
      P = MR(1)
      Q = MR(3)
      PS = MR(4)
      QS = MR(6)
      IS = MR(7)
      IFAULT = 0
      IF (P.LT.0) IFAULT = 1
      IF (Q.LT.0) IFAULT = 1
      IF (PS.LT.0) IFAULT = 1
      IF (QS.LT.0) IFAULT = 1
      IF (IS.LT.0) IFAULT = 1
      IF ((IS.EQ.0) .AND. (PS.GT.0 .OR. QS.GT.0)) IFAULT = 1
      IF (NPAR.EQ.0 .OR. NPAR.NE.(P+Q+PS+QS)) IFAULT = 1
      IF (M.LE.NPAR) IFAULT = 1
      IF (M.GE.N) IFAULT = 1
      IF (IM.LT.M) IFAULT = 1
      IF (IFAULT.NE.0) THEN
         WRITE (P01REC,FMT=99999) N, P, Q, PS, QS, IS, M, IM, NPAR
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,4,P01REC)
         RETURN
      END IF
C
      SHOW = .FALSE.
      IF (ISHOW.GT.0) SHOW = .TRUE.
      ISEA = IS
      L2 = MAX(P,Q,PS,QS)
      IF (LMAX.LT.L2) THEN
         IFAULT = 1
         WRITE (P01REC,FMT=99998) LMAX, L2
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
         RETURN
      END IF
C
      SEASON = .TRUE.
      IF (IS.EQ.0) THEN
         SEASON = .FALSE.
         ISEA = 1
      END IF
C
      I2 = NPAR*(M+NPAR+1) + L2*ISEA + M
      IF (LWORK.LT.I2) THEN
         IFAULT = 1
         WRITE (P01REC,FMT=99997) LWORK, I2
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
         RETURN
      END IF
C
C     test for stationarity and invertibility
C
      K = 1
      IFAULT = 2
      IF (P.GT.0) THEN
         LW2 = P*P + 1
         LW3 = LW2 + P
         LIMIT = 1.0D0
         CALL G05HDZ(P,K,PAR,WORK,WORK(LW2),WORK(LW3),INTGR,P,LIMIT,
     *               BIGEST,STAT,IERR)
         IF (IERR.EQ.1) IFAULT = 4
         IF (( .NOT. STAT) .OR. (IERR.EQ.1)) THEN
            IF (IFAULT.EQ.2) THEN
               WRITE (P01REC,FMT=99996)
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
            ELSE
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,0,P01REC)
            END IF
            RETURN
         END IF
      END IF
C
      IF (Q.GT.0) THEN
         LW2 = Q*Q + 1
         LW3 = LW2 + Q
         LIMIT = 1.0D0
         CALL G05HDZ(Q,K,PAR(P+1),WORK,WORK(LW2),WORK(LW3),INTGR,Q,
     *               LIMIT,BIGEST,STAT,IERR)
         IF (IERR.EQ.1) IFAULT = 4
         IF (( .NOT. STAT) .OR. (IERR.EQ.1)) THEN
            IF (IFAULT.EQ.2) THEN
               WRITE (P01REC,FMT=99995)
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
            ELSE
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,0,P01REC)
            END IF
            RETURN
         END IF
      END IF
C
      IF (PS.GT.0) THEN
         LW2 = PS*PS + 1
         LW3 = LW2 + PS
         LIMIT = 1.0D0
         CALL G05HDZ(PS,K,PAR(P+Q+1),WORK,WORK(LW2),WORK(LW3),INTGR,PS,
     *               LIMIT,BIGEST,STAT,IERR)
         IF (IERR.EQ.1) IFAULT = 4
         IF (( .NOT. STAT) .OR. (IERR.EQ.1)) THEN
            IF (IFAULT.EQ.2) THEN
               WRITE (P01REC,FMT=99994)
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
            ELSE
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,0,P01REC)
            END IF
            RETURN
         END IF
      END IF
C
      IF (QS.GT.0) THEN
         LW2 = QS*QS + 1
         LW3 = LW2 + QS
         LIMIT = 1.0D0
         CALL G05HDZ(QS,K,PAR(P+Q+PS+1),WORK,WORK(LW2),WORK(LW3),INTGR,
     *               QS,LIMIT,BIGEST,STAT,IERR)
         IF (IERR.EQ.1) IFAULT = 4
         IF (( .NOT. STAT) .OR. (IERR.EQ.1)) THEN
            IF (IFAULT.EQ.2) THEN
               WRITE (P01REC,FMT=99993)
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
            ELSE
               IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,0,P01REC)
            END IF
            RETURN
         END IF
      END IF
      IFAULT = 0
C
C     calculate residual autocorrelation function
C
      IFAILA = 1
      CALL G13ABF(V,N,M,XM,XV,C,QSTAT,IFAILA)
      IF (IFAILA.GT.0) THEN
         IFAULT = 3
         DO 20 I = 1, M
            C(I) = 0.0D0
   20    CONTINUE
      END IF
C
C     calculate the Box - Ljung portmanteau statistic
C
      SUM2 = 0.0D0
      DO 40 I = 1, M
         SUM2 = SUM2 + C(I)*C(I)/DBLE(N-I)
   40 CONTINUE
C
      SUM2 = SUM2*DBLE(N)*DBLE(N+2)
      IDF = M - P - Q - PS - QS
      IFAILA = 1
      SIGLEV = G01ECF('U',SUM2,DBLE(IDF),IFAILA)
      IF (IFAULT.EQ.3) GO TO 220
C
C                                       -1
C     call G13ASX to calculate X   (X'X)   X'
C
      LW2 = M*NPAR + NPAR*(NPAR+1) + 1
C
      CALL G13ASX(PAR,P,Q,PS,QS,WORK,WORK(M*NPAR+1),ACFVAR,NPAR,M,IM,
     *            ISEA,WORK(LW2),N,L2,IERR)
C
      IF (IERR.GT.0) THEN
C
C        set standard errors to 1 / sqrt(n)
C
         SUM = 1.0D0/SQRT(DBLE(N))
         DO 80 J = 1, M
            DO 60 I = 1, M
               ACFVAR(I,J) = 0.0D0
   60       CONTINUE
            ACFVAR(J,J) = SUM
   80    CONTINUE
         IF (IFAULT.EQ.0) IFAULT = 5
C
      ELSE
C
C        convert ACFVAR to a standard error / correlation matrix
C
         I2 = 0
         DO 100 I = 1, M
            IF (ACFVAR(I,I).GT.0.0D0) THEN
               ACFVAR(I,I) = SQRT(ACFVAR(I,I))
               I2 = I2 + 1
            END IF
  100    CONTINUE
C
         IF (I2.EQ.M) THEN
            DO 140 I = 2, M
               SUM = ACFVAR(I,I)
               DO 120 J = 1, I - 1
                  ACFVAR(I,J) = ACFVAR(I,J)/(SUM*ACFVAR(J,J))
                  ACFVAR(J,I) = ACFVAR(I,J)
  120          CONTINUE
  140       CONTINUE
         ELSE
            SUM = 1.0D0/SQRT(DBLE(N))
            DO 180 J = 1, M
               DO 160 I = 1, M
                  ACFVAR(I,J) = 0.0D0
  160          CONTINUE
               ACFVAR(J,J) = SUM
  180       CONTINUE
            IF (IFAULT.EQ.0) IFAULT = 6
         END IF
      END IF
C
      IF (SHOW) THEN
C
         CALL X04ABF(0,LP)
         WRITE (REC,FMT=99992)
         CALL X04BAY(LP,4,REC)
         L2 = M/7
         DO 200 L = 1, L2
            I = (L-1)*7
            WRITE (REC,FMT=99991) (I+J,J=1,7)
            CALL X04BAF(LP,REC(1))
            WRITE (REC,FMT=99990) (C(I+J),J=1,7)
            CALL X04BAF(LP,REC(1))
            WRITE (REC,FMT=99989) (ACFVAR(I+J,I+J),J=1,7)
            CALL X04BAF(LP,REC(1))
            WRITE (REC,FMT=99988)
            CALL X04BAF(LP,REC(1))
  200    CONTINUE
         L3 = 7*(M/7)
         IF (M.GT.L3) THEN
            WRITE (REC,FMT=99991) (I,I=L3+1,M)
            CALL X04BAF(LP,REC(1))
            WRITE (REC,FMT=99990) (C(I),I=L3+1,M)
            CALL X04BAF(LP,REC(1))
            WRITE (REC,FMT=99989) (ACFVAR(I,I),I=L3+1,M)
            CALL X04BAF(LP,REC(1))
            WRITE (REC,FMT=99988)
            CALL X04BAF(LP,REC(1))
         END IF
C
         WRITE (REC,FMT=99987) SUM2, SIGLEV, IDF
         CALL X04BAY(LP,4,REC)
      END IF
C
  220 IF (IFAULT.EQ.3) THEN
         WRITE (P01REC,FMT=99986)
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
      ELSE IF (IFAULT.EQ.5) THEN
         WRITE (P01REC,FMT=99985)
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,2,P01REC)
      ELSE IF (IFAULT.EQ.6) THEN
         WRITE (P01REC,FMT=99984)
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,2,P01REC)
      ELSE
         IFAIL = 0
      END IF
C
      IF (SHOW) THEN
         WRITE (REC,FMT=99983) IFAIL
         CALL X04BAY(LP,2,REC)
      END IF
      RETURN
C
99999 FORMAT ('  ** ON ENTRY, ONE OR MORE OF THE FOLLOWING PARAMETER V',
     *       'ALUES IS ILLEGAL',/'  N    = ',I16,' MR(1) = ',I16,' MR(',
     *       '3) = ',I16,/' MR(4) = ',I16,' MR(6) = ',I16,' MR(7) = ',
     *       I16,/'  M    = ',I16,'  IRCM  = ',I16,'  NPAR = ',I16)
99998 FORMAT ('  ** ON ENTRY, LIW ( ',I16,' ) MUST BE AT LEAST ',I16)
99997 FORMAT ('  ** ON ENTRY, LWORK ( ',I16,' ) MUST BE AT LEAST ',I16)
99996 FORMAT ('  ** ON ENTRY, THE NONSEASONAL AUTOREGRESSIVE OPERATOR ',
     *       'IS NONSTATIONARY')
99995 FORMAT ('  ** ON ENTRY, THE NONSEASONAL MOVING AVERAGE OPERATOR ',
     *       'IS NON-INVERTIBLE')
99994 FORMAT ('  ** ON ENTRY, THE SEASONAL AUTOREGRESSIVE OPERATOR IS ',
     *       'NONSTATIONARY')
99993 FORMAT ('  ** ON ENTRY, THE SEASONAL MOVING AVERAGE OPERATOR IS ',
     *       'NON-INVERTIBLE')
99992 FORMAT (/' RESIDUAL AUTOCORRELATION FUNCTION',/' ---------------',
     *       '------------------',/)
99991 FORMAT (' LAG  K  ',I5,6I7)
99990 FORMAT (' R(K)    ',7F7.3)
99989 FORMAT (' ST.ERROR',7F7.3)
99988 FORMAT (' ------------------------------------------------------',
     *       '---')
99987 FORMAT (/' BOX - LJUNG PORTMANTEAU STATISTIC = ',F10.3,/16X,'SIG',
     *       'NIFICANCE LEVEL = ',F10.3,/' (BASED ON ',I3,' DEGREES OF',
     *       ' FREEDOM)')
99986 FORMAT ('  ** ON ENTRY, THE ELEMENTS OF V ARE NEARLY IDENTICAL G',
     *       'IVING NEAR-ZERO VARIANCE')
99985 FORMAT ('  ** ON ENTRY, ONE OR MORE OF THE AR OPERATORS HAS A FA',
     *       'CTOR IN COMMON WITH ONE',/'     OR MORE OF THE MA OPERAT',
     *       'ORS')
99984 FORMAT ('  ** THE MATRIX RCM COULD NOT BE COMPUTED BECAUSE',/'  ',
     *       '   ONE OF ITS DIAGONAL ELEMENTS WAS FOUND TO BE NON-POSI',
     *       'TIVE')
99983 FORMAT (/' VALUE OF IFAIL PARAMETER ON EXIT FROM G13ASF = ',I3)
      END
