      SUBROUTINE G13DSF(K,N,V,IK,P,Q,M,PAR,PARHLD,QQ,ISHOW,R0,C,ACFVAR,
     *                  IM,CHI,IDF,SIGLEV,IW,LIW,WORK,LWORK,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DSF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHI, SIGLEV
      INTEGER           IDF, IFAIL, IK, IM, ISHOW, K, LIW, LWORK, M, N,
     *                  P, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  ACFVAR(IM,M*K*K), C(IK,IK,M), PAR((P+Q)*K*K),
     *                  QQ(IK,K), R0(IK,K), V(IK,N), WORK(LWORK)
      INTEGER           IW(LIW)
      LOGICAL           PARHLD((P+Q)*K*K)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGEST, LIMIT, SUM2
      INTEGER           I, I2, IERR, IFAILA, IFAULT, II, ITW, J, K2, KP,
     *                  KQ, L, L1, L2, L3, LP, LW1, LW2, LW3, LW4, LW5,
     *                  LW6, LW7, LW8, MK2, MM, NPAR
      LOGICAL           NOPARS, SHOW, STAT, VEQUAL
C     .. Local Arrays ..
      CHARACTER         ST(80)
      CHARACTER*80      P01REC(4), REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03ABF, G05HDZ, G13DMF, G13DST, G13DSU, G13DSV,
     *                  G13DSW, G13DSZ, X04ABF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MOD, DBLE, SQRT
C     .. Executable Statements ..
C
C     ******************************************************************
C
C     This routine calculates the variance-covariance matrix of the
C     residual cross-correlation matrices in a vector ARMA model.
C     It also returns a portmanteau statistic along with its'
C     significance level. A user would normally summon this routine
C     after calling G13DCF to fit a vector ARMA model.
C
C     The computational procedure is described in a paper by Li, W.K.
C     and McLeod, A.I.
C         'Distribution of the residual autocorrelations in multivariate
C          ARMA time series models'
C          Journal of the Royal Statistical Society, Series B, Vol 43,
C          pages 231 - 239, 1981.
C
C     ******************************************************************
C
C     first test for errors in the input arguments
C
      SHOW = .FALSE.
      IF (ISHOW.NE.0) SHOW = .TRUE.
      IFAULT = 0
      IF (K.LT.1) IFAULT = 1
      IF (IK.LT.K) IFAULT = 1
      IF (P.LT.0) IFAULT = 1
      IF (Q.LT.0) IFAULT = 1
      NPAR = (P+Q)*K*K
      IF (NPAR.EQ.0) IFAULT = 1
      IF (M.LE.P+Q) IFAULT = 1
      IF (M.GE.N) IFAULT = 1
      IF (IM.LT.M*K*K) IFAULT = 1
      IF (IFAULT.NE.0) THEN
         WRITE (P01REC,FMT=99998) K, N, IK, P, Q, M, IM
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,4,P01REC)
         RETURN
      END IF
C
C     set NOPARS depending on whether all elements of PARHLD are
C     equal to .TRUE.
C
      NOPARS = .FALSE.
      DO 20 I = 1, NPAR
         IF ( .NOT. PARHLD(I)) GO TO 40
   20 CONTINUE
      NOPARS = .TRUE.
C
C     test whether the workspace arrays are big enough
C
   40 J = K*MAX(P,Q)
      IF (LIW.LT.J) THEN
         IFAULT = 1
         WRITE (P01REC,FMT=99997) LIW, J
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
         RETURN
      END IF
C
      MK2 = M*K*K
      J = 2*K + MK2*(NPAR+MK2+1) + 3*K*K + (NPAR+1)*NPAR
C
      IF (LWORK.LT.J) THEN
         IFAULT = 1
         WRITE (P01REC,FMT=99996) LWORK, J
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
         RETURN
      END IF
C
C     test whether QQ is positive-definite
C
C     first set the upper triangle of QQ to the lower triangle of QQ
C
      DO 80 I = 1, K
         DO 60 J = 1, I
            QQ(J,I) = QQ(I,J)
   60    CONTINUE
   80 CONTINUE
C
      IFAILA = 1
      CALL F03ABF(QQ,IK,K,SUM2,WORK,IFAILA)
C
C     reset the strict lower triangle of QQ (which F03ABF has corrupted)
C     to its former value
C
      DO 120 I = 1, K
         DO 100 J = 1, I - 1
            QQ(I,J) = QQ(J,I)
  100    CONTINUE
  120 CONTINUE
C
      IF (IFAILA.GT.0) THEN
         IFAULT = 2
         WRITE (P01REC,FMT=99995)
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
         RETURN
      END IF
C
      K2 = K*K
C
C     test for stationarity of the AR operator
C
      IFAULT = 2
      IF (P.GT.0) THEN
         KP = P*K
         LW1 = 1
         LW2 = LW1 + KP*KP
         LW3 = LW2 + KP
         LIMIT = 1.0D0
         CALL G05HDZ(P,K,PAR,WORK(LW1),WORK(LW2),WORK(LW3),IW,KP,LIMIT,
     *               BIGEST,STAT,IERR)
         IF (IERR.EQ.1) IFAULT = 4
         IF (( .NOT. STAT) .OR. IERR.EQ.1) THEN
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
C     test for invertibility of the MA operator
C
      IF (Q.GT.0) THEN
         KQ = K*Q
         LW1 = 1
         LW2 = LW1 + KQ*KQ
         LW3 = LW2 + KQ
         LIMIT = 1.0D0
         CALL G05HDZ(Q,K,PAR(P*K*K+1),WORK(LW1),WORK(LW2),WORK(LW3),IW,
     *               KQ,LIMIT,BIGEST,STAT,IERR)
         IF (IERR.EQ.1) IFAULT = 4
         IF (( .NOT. STAT) .OR. IERR.EQ.1) THEN
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
C     calculate the residual cross correlation matrices
C
      VEQUAL = .FALSE.
      IFAILA = 1
      CALL G13DMF('R-correlation',K,N,M,V,IK,WORK,R0,C,IFAILA)
      IF (IFAILA.GT.0) IFAULT = 3
C
C     set the diagonal elements of R0 to one and store in first k * k
C     elements of WORK array
C
      DO 160 I = 1, K
         DO 140 J = 1, K
            WORK((J-1)*K+I) = R0(I,J)
  140    CONTINUE
         WORK((I-1)*K+I) = 1.0D0
  160 CONTINUE
C
C     calculate Li and McLeod's modified portmanteau statistic
C
      LW1 = K2 + 1
      LW2 = LW1 + (K2+1)*K2
C
      CALL G13DSZ(IK,K,K2,M,N,C,P,Q,PARHLD,NPAR,WORK,WORK(LW1),WORK(LW2)
     *            ,CHI,IDF,SIGLEV,IFAILA,IFAULT)
      IF (IFAILA.GT.0) THEN
         IFAULT = 3
         VEQUAL = .TRUE.
      END IF
      IF (IFAULT.EQ.3) THEN
         CHI = 0.0D0
         SIGLEV = 1.0D0
      END IF
C
C     call a subroutine to calculate the elements of ACFVAR
C     but first partition the workspace array WORK
C
      MK2 = M*K2
      LW1 = K2 + 1
      LW2 = LW1 + NPAR*MK2
      LW3 = LW2 + K
      LW4 = LW3 + MK2
      LW5 = LW4 + K2
      LW6 = LW5 + K
      LW7 = LW6 + (NPAR+1)*NPAR
      LW8 = LW7 + K2
C
      CALL G13DSW(K,P,Q,M,QQ,IK,ACFVAR,IM,WORK,WORK(LW1),WORK(LW2),
     *            WORK(LW3),WORK(LW4),N,PAR,NPAR,PARHLD,WORK(LW5),
     *            WORK(LW6),WORK(LW7),WORK(LW8),NOPARS,IERR)
C
      IF (IERR.GT.0) THEN
         SUM2 = 1.0D0/SQRT(DBLE(N))
         DO 200 J = 1, MK2
            DO 180 I = 1, MK2
               ACFVAR(I,J) = 0.0D0
  180       CONTINUE
            ACFVAR(J,J) = SUM2
  200    CONTINUE
         IF (IFAULT.EQ.0) IFAULT = 5
      ELSE
C
C        convert ACFVAR to a standard error - correlation matrix
C
         I2 = 0
         DO 220 I = 1, MK2
            IF (ACFVAR(I,I).GT.0.0D0) THEN
               ACFVAR(I,I) = SQRT(ACFVAR(I,I))
               I2 = I2 + 1
            END IF
  220    CONTINUE
         IF (I2.EQ.MK2) THEN
            DO 260 I = 2, MK2
               SUM2 = ACFVAR(I,I)
               DO 240 J = 1, I - 1
                  ACFVAR(I,J) = ACFVAR(I,J)/(SUM2*ACFVAR(J,J))
                  ACFVAR(J,I) = ACFVAR(I,J)
  240          CONTINUE
  260       CONTINUE
         ELSE
            SUM2 = 1.0D0/SQRT(DBLE(N))
            DO 300 J = 1, MK2
               DO 280 I = 1, MK2
                  ACFVAR(I,J) = 0.0D0
  280          CONTINUE
               ACFVAR(J,J) = SUM2
  300       CONTINUE
            IF (IFAULT.EQ.0) IFAULT = 6
         END IF
      END IF
C
      IF (SHOW) THEN
         CALL X04ABF(0,LP)
C
         IF (K.EQ.1) WRITE (REC,FMT=99999)
         IF (K.GT.1) WRITE (REC,FMT=99992)
         CALL X04BAY(LP,3,REC)
         DO 340 L = 1, M
            WRITE (REC,FMT=99991) L, (C(1,J,L),J=1,K)
            CALL X04BAY(LP,2,REC)
            L2 = (L-1)*K2 + 1
            WRITE (REC,FMT=99989) (ACFVAR(L2+(J-1)*K,L2+(J-1)*K),J=1,K)
            CALL X04BAF(LP,REC(1))
            DO 320 I = 2, K
               WRITE (REC,FMT=99990) (C(I,J,L),J=1,K)
               CALL X04BAF(LP,REC(1))
               L2 = (L-1)*K2 + I
               WRITE (REC,FMT=99989) (ACFVAR(L2+(J-1)*K,L2+(J-1)*K),J=1,
     *           K)
               CALL X04BAF(LP,REC(1))
  320       CONTINUE
  340    CONTINUE
C
         IF (K.LE.6) THEN
            WRITE (REC,FMT=99988)
            CALL X04BAY(LP,3,REC)
            MM = (79/K) - 3
            DO 380 II = 1, M/MM
               ITW = (MM+3)*K + 1
               L1 = (II-1)*MM + 1
               L3 = L1 + MM - 1
               WRITE (REC,FMT=99987) L1, L3
               CALL X04BAY(LP,3,REC)
C
C              print out a row of *'s
C
               CALL G13DST(ST,ITW,LP)
C
               DO 360 I = 1, K
C
C                 now print out significant cross-correlations
C
                  I2 = I
                  CALL G13DSU(ST,ITW,LP,MM,K)
                  CALL G13DSV(K,I2,ST,L1,L3,ITW,MM,ACFVAR,IM,IK,C,LP,M)
                  CALL G13DSU(ST,ITW,LP,MM,K)
                  CALL G13DST(ST,ITW,LP)
C
  360          CONTINUE
C
  380       CONTINUE
C
            IF (MOD(M,MM).GT.0) THEN
               ITW = M - MM*(M/MM)
               L1 = M - ITW + 1
               L3 = M
               MM = ITW
               ITW = (MM+3)*K + 1
               WRITE (REC,FMT=99987) L1, L3
               CALL X04BAY(LP,3,REC)
C
C              print out a row of *'s
C
               CALL G13DST(ST,ITW,LP)
C
               DO 400 I = 1, K
C
C                 now print out significant cross-correlations
C
                  I2 = I
                  CALL G13DSU(ST,ITW,LP,MM,K)
                  CALL G13DSV(K,I2,ST,L1,L3,ITW,MM,ACFVAR,IM,IK,C,LP,M)
                  CALL G13DSU(ST,ITW,LP,MM,K)
                  CALL G13DST(ST,ITW,LP)
C
  400          CONTINUE
C
            END IF
         END IF
C
         WRITE (REC,FMT=99986) CHI, SIGLEV, IDF
         CALL X04BAY(LP,4,REC)
C
      END IF
C
      IF (IFAULT.EQ.3) THEN
         IF (VEQUAL) THEN
            WRITE (P01REC,FMT=99985)
            IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
         ELSE
            WRITE (P01REC,FMT=99984)
            IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,2,P01REC)
         END IF
      ELSE IF (IFAULT.EQ.5) THEN
         WRITE (P01REC,FMT=99983)
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,1,P01REC)
      ELSE IF (IFAULT.EQ.6) THEN
         WRITE (P01REC,FMT=99982)
         IFAIL = P01ABF(IFAIL,IFAULT,SRNAME,2,P01REC)
      ELSE
         IFAIL = 0
      END IF
C
      IF (SHOW) THEN
         WRITE (REC,FMT=99981) IFAIL
         CALL X04BAY(LP,2,REC)
      END IF
      RETURN
C
99999 FORMAT (/' RESIDUAL AUTOCORRELATION FUNCTION',/' ---------------',
     *       '------------------')
99998 FORMAT ('  ** ON ENTRY, ONE OR MORE OF THE FOLLOWING PARAMETER V',
     *       'ALUES IS ILLEGAL',/' K   = ',I16,'   N  = ',I16,
     *       '   IK = ',I16,/' IP  = ',I16,'   IQ = ',I16,'   M  = ',
     *       I16,/' IRCM = ',I16)
99997 FORMAT ('  ** ON ENTRY, LIW ( ',I16,' ) MUST BE AT LEAST ',I16)
99996 FORMAT ('  ** ON ENTRY, LWORK (',I16,' ) MUST BE AT LEAST ',I16)
99995 FORMAT ('  ** ON ENTRY, THE COVARIANCE MATRIX QQ IS NOT POSITIVE',
     *       '-DEFINITE')
99994 FORMAT ('  ** ON ENTRY, THE AR PARAMETER ESTIMATES ARE OUTSIDE T',
     *       'HE STATIONARITY REGION')
99993 FORMAT ('  ** ON ENTRY, THE MA PARAMETER ESTIMATES ARE OUTSIDE T',
     *       'HE INVERTIBILITY REGION')
99992 FORMAT (/' RESIDUAL CROSS-CORRELATION MATRICES',/' -------------',
     *       '----------------------')
99991 FORMAT (/' LAG   ',I3,'           : ',7F8.3)
99990 FORMAT (23X,7F8.3)
99989 FORMAT (24X,7(' (',F5.3,')',:))
99988 FORMAT (/' SUMMARY TABLE',/' -------------')
99987 FORMAT (/' LAGS ',I3,' - ',I3,/)
99986 FORMAT (/' LI-MCLEOD PORTMANTEAU STATISTIC = ',F10.3,/14X,'SIGNI',
     *       'FICANCE LEVEL = ',F10.3,/' (BASED ON ',I3,' DEGREES OF F',
     *       'REEDOM)')
99985 FORMAT ('  ** ON ENTRY, AT LEAST TWO OF THE RESIDUAL SERIES ARE ',
     *       'IDENTICAL')
99984 FORMAT ('  ** ON ENTRY, AT LEAST ONE OF THE RESIDUAL SERIES IN T',
     *       'HE ARRAY V',/'     HAS NEAR-ZERO VARIANCE')
99983 FORMAT ('  ** ON ENTRY, THE AR OPERATOR HAS A FACTOR IN COMMON W',
     *       'ITH THE MA OPERATOR')
99982 FORMAT ('  ** THE MATRIX RCM COULD NOT BE COMPUTED BECAUSE',/'  ',
     *       '   ONE OF ITS DIAGONAL ELEMENTS WAS FOUND TO BE NON-POSI',
     *       'TIVE')
99981 FORMAT (/' VALUE OF IFAIL PARAMETER ON EXIT FROM G13DSF = ',I2)
      END
