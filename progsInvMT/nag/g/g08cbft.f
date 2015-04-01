      SUBROUTINE G08CBF(N,X,DIST,PAR,ESTIMA,NTYPE,D,Z,P,SX,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08CBF')
      DOUBLE PRECISION  BIG
      PARAMETER         (BIG=1.0D06)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D, P, Z
      INTEGER           IFAIL, N, NTYPE
      CHARACTER*1       ESTIMA
      CHARACTER*(*)     DIST
C     .. Array Arguments ..
      DOUBLE PRECISION  PAR(2), SX(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DN, DP, DT, EN, FF, FN, FO, PDF, PEQK, PEXP,
     *                  PGTK, Q, SR, SUM, SUMSQ, SVAR, TOL, XBAR, XMAX,
     *                  XMIN, XNORM, Y
      INTEGER           I, IERROR, IF2, J, KNEW, KPREV, NBIN, NREC
      LOGICAL           BETA, BINOM, DISCRE, EXPON, FIRST4, GAMMA, LEST,
     *                  NORMAL, ONEPAR, POISON, TWOPAR, UNIFOR
C     .. Local Arrays ..
      CHARACTER*80      P01REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  G08CBZ, S15ABF, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          G08CBZ, S15ABF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01BJF, G01BKF, G01EEF, M01CAF, S14BAF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, INT, LEN, LOG, MAX, DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 0
      IERROR = 0
      LEST = ESTIMA .EQ. 'E' .OR. ESTIMA .EQ. 'e'
      BETA = .FALSE.
      BINOM = .FALSE.
      UNIFOR = DIST(1:1) .EQ. 'U' .OR. DIST(1:1) .EQ. 'u'
      NORMAL = DIST(1:1) .EQ. 'N' .OR. DIST(1:1) .EQ. 'n'
      GAMMA = DIST(1:1) .EQ. 'G' .OR. DIST(1:1) .EQ. 'g'
      IF (LEN(DIST).GE.2) THEN
         BETA = DIST(1:2) .EQ. 'BE' .OR. DIST(1:2)
     *          .EQ. 'Be' .OR. DIST(1:2) .EQ. 'bE' .OR. DIST(1:2)
     *          .EQ. 'be'
         BINOM = DIST(1:2) .EQ. 'BI' .OR. DIST(1:2)
     *           .EQ. 'Bi' .OR. DIST(1:2) .EQ. 'bI' .OR. DIST(1:2)
     *           .EQ. 'bi'
      END IF
      EXPON = DIST(1:1) .EQ. 'E' .OR. DIST(1:1) .EQ. 'e'
      POISON = DIST(1:1) .EQ. 'P' .OR. DIST(1:1) .EQ. 'p'
C
C     Set up some other useful logicals
C
      TWOPAR = UNIFOR .OR. NORMAL .OR. GAMMA .OR. BETA .OR. BINOM
      ONEPAR = EXPON .OR. POISON
      DISCRE = BINOM .OR. POISON
      FIRST4 = UNIFOR .OR. NORMAL .OR. GAMMA .OR. BETA
C
C     Start error checks
C
      IF (N.LT.3) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF ( .NOT. UNIFOR .AND. .NOT. NORMAL .AND. .NOT. GAMMA .AND.
     *         .NOT. BETA .AND. .NOT. BINOM .AND. .NOT. EXPON .AND.
     *         .NOT. POISON) THEN
         IERROR = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) DIST(1:1)
      ELSE IF (NTYPE.LE.0 .OR. NTYPE.GT.3) THEN
         IERROR = 3
         NREC = 1
         WRITE (P01REC,FMT=99997) NTYPE
      ELSE IF ( .NOT. LEST .AND. ESTIMA.NE.'S' .AND. ESTIMA.NE.'s') THEN
         IERROR = 4
         NREC = 1
         WRITE (P01REC,FMT=99996) ESTIMA
      ELSE IF ( .NOT. LEST .AND. UNIFOR .AND. PAR(1).GE.PAR(2)) THEN
         IERROR = 5
         NREC = 3
         WRITE (P01REC,FMT=99995) PAR(1), PAR(2)
      ELSE IF ( .NOT. LEST .AND. NORMAL .AND. PAR(2).LE.0.0D0) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99994) PAR(2)
      ELSE IF ( .NOT. LEST .AND. GAMMA .AND. (PAR(1)
     *         .LE.0.0D0 .OR. PAR(2).LE.0.0D0)) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99993) PAR(1), PAR(2)
      ELSE IF ( .NOT. LEST .AND. BETA .AND. (PAR(1).LE.0.0D0 .OR. PAR(1)
     *         .GT.BIG .OR. PAR(2).LE.0.0D0 .OR. PAR(2).GT.BIG)) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99992) PAR(1), PAR(2)
      ELSE IF (BINOM .AND. (PAR(1).LT.1.0D0 .OR. (PAR(1)*X02AJF())
     *         .GT.1.0D0)) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99991) PAR(1)
      ELSE IF ( .NOT. LEST .AND. BINOM .AND. (PAR(2)
     *         .LE.0.0D0 .OR. PAR(2).GE.1.0D0)) THEN
         IERROR = 5
         NREC = 3
         WRITE (P01REC,FMT=99990) PAR(2)
      ELSE IF ( .NOT. LEST .AND. EXPON .AND. PAR(1).LE.0.0D0) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99989) PAR(1)
      ELSE IF ( .NOT. LEST .AND. POISON .AND. (PAR(1)
     *         .LE.0.0D0 .OR. PAR(1).GT.BIG)) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99988) PAR(1)
      ELSE
         EN = DBLE(N)
         DO 20 I = 1, N
            SX(I) = X(I)
   20    CONTINUE
         CALL M01CAF(SX,1,N,'A',IF2)
         XMIN = SX(1)
         XMAX = SX(N)
C
C        Check data against the support of the chosen distribution
C
         IF (UNIFOR .AND. .NOT. LEST .AND. (XMIN.LT.PAR(1)
     *       .OR. XMAX.GT.PAR(2))) THEN
            IERROR = 6
            NREC = 3
            WRITE (P01REC,FMT=99987) PAR(1), PAR(2)
         ELSE IF (GAMMA .AND. XMIN.LT.0.0D0) THEN
            IERROR = 6
            NREC = 2
            WRITE (P01REC,FMT=99986)
         ELSE IF (BETA .AND. (XMIN.LT.0.0D0 .OR. XMAX.GT.1.0D0)) THEN
            IERROR = 6
            NREC = 3
            WRITE (P01REC,FMT=99985)
         ELSE IF (BINOM .AND. (XMIN.LT.0.0D0 .OR. XMAX.GT.PAR(1))) THEN
            IERROR = 6
            NREC = 4
            WRITE (P01REC,FMT=99984) PAR(1)
         ELSE IF (BINOM .AND. LEST .AND. XMAX.LE.0.0D0) THEN
            IERROR = 6
            NREC = 2
            WRITE (P01REC,FMT=99978)
         ELSE IF (EXPON .AND. XMIN.LT.0.0D0) THEN
            IERROR = 6
            NREC = 2
            WRITE (P01REC,FMT=99983)
         ELSE IF (EXPON .AND. LEST .AND. XMAX.LE.0.0D0) THEN
            IERROR = 6
            NREC = 2
            WRITE (P01REC,FMT=99977)
         ELSE IF (POISON .AND. XMIN.LT.0.0D0) THEN
            IERROR = 6
            NREC = 2
            WRITE (P01REC,FMT=99982)
         ELSE IF (POISON .AND. LEST .AND. XMAX.LE.0.0D0) THEN
            IERROR = 6
            NREC = 2
            WRITE (P01REC,FMT=99976)
         ELSE IF (FIRST4 .AND. (XMIN.GE.XMAX) .AND. LEST) THEN
            IERROR = 7
            NREC = 2
            WRITE (P01REC,FMT=99981)
         ELSE
C
C              Estimate parameters if necessary
C
            IF (LEST) THEN
               IF (UNIFOR) THEN
                  PAR(1) = XMIN
                  PAR(2) = XMAX
               ELSE
                  SUM = 0.0D0
                  DO 40 I = 1, N
                     SUM = SUM + X(I)
   40             CONTINUE
                  XBAR = SUM/EN
                  IF (NORMAL .OR. POISON) PAR(1) = XBAR
                  IF (BINOM) PAR(2) = XBAR/PAR(1)
                  IF (EXPON) PAR(1) = 1.0D0/XBAR
                  IF (FIRST4) THEN
                     SUMSQ = 0.0D0
                     DO 60 I = 1, N
                        SUMSQ = SUMSQ + (X(I)-XBAR)*(X(I)-XBAR)
   60                CONTINUE
                     SVAR = SUMSQ/(EN-1.0D0)
                     IF (NORMAL) THEN
                        PAR(2) = SVAR
                     ELSE IF (GAMMA) THEN
                        PAR(1) = SVAR/XBAR
                        PAR(2) = (XBAR**2)/SVAR
                     ELSE IF (BETA) THEN
                        PAR(1) = XBAR*((XBAR*(1.0D0-XBAR)/SVAR)-1.0D0)
                        PAR(2) = (1.0D0-XBAR)*((XBAR*(1.0D0-XBAR)/SVAR)
     *                           -1.0D0)
                     END IF
                  END IF
               END IF
            END IF
            IF (BINOM) THEN
               IF ((PAR(1)*PAR(2)*(1.0D0-PAR(2))).GT.BIG) THEN
                  IERROR = 8
                  NREC = 2
                  WRITE (P01REC,FMT=99980) PAR(1), PAR(2)
                  GO TO 120
               END IF
               NBIN = PAR(1)
            END IF
            IF (GAMMA .OR. BETA) TOL = 100.0D0*X02AJF()
            IF (EXPON) SR = -LOG(X02AMF())
            FO = 0.0D0
            KPREV = 0
            D = 0.0D0
            DP = 0.0D0
            DN = 0.0D0
            DO 100 J = 1, N
               FN = DBLE(J)/EN
               IF (UNIFOR) THEN
                  FF = (SX(J)-PAR(1))/(PAR(2)-PAR(1))
                  IF (FF.LT.0.0D0) FF = 0.0D0
                  IF (FF.GT.1.0D0) FF = 1.0D0
               ELSE IF (NORMAL) THEN
                  XNORM = (SX(J)-PAR(1))/SQRT(PAR(2))
                  IF2 = 0
                  FF = S15ABF(XNORM,IF2)
               ELSE IF (GAMMA) THEN
                  IF2 = 1
                  Y = SX(J)/PAR(2)
                  CALL S14BAF(PAR(1),Y,TOL,FF,Q,IF2)
                  IF (IF2.EQ.3) THEN
                     IERROR = 9
                     NREC = 3
                     WRITE (P01REC,FMT=99979)
                     GO TO 120
                  END IF
               ELSE IF (BETA) THEN
                  IF2 = 1
                  CALL G01EEF(SX(J),PAR(1),PAR(2),TOL,FF,Q,PDF,IF2)
               ELSE IF (BINOM) THEN
                  KNEW = INT(SX(J))
                  IF (KNEW.EQ.KPREV .AND. J.LT.N) THEN
                     GO TO 80
                  ELSE
                     IF2 = 0
                     CALL G01BJF(NBIN,PAR(2),KPREV,FF,PGTK,PEQK,IF2)
                     KPREV = KNEW
                  END IF
               ELSE IF (EXPON) THEN
                  PEXP = PAR(1)*SX(J)
                  IF (PEXP.LT.SR) THEN
                     FF = 1.0D0 - EXP(-PEXP)
                  ELSE
                     FF = 1.0D0
                  END IF
               ELSE IF (POISON) THEN
                  KNEW = INT(SX(J))
                  IF (KNEW.EQ.KPREV) THEN
                     GO TO 80
                  ELSE
                     IF2 = 0
                     CALL G01BKF(PAR(1),KPREV,FF,PGTK,PEQK,IF2)
                     KPREV = KNEW
                  END IF
               END IF
               IF (DISCRE) THEN
                  DT = FO - FF
               ELSE
                  DT = FN - FF
               END IF
               IF (DT.GT.DP) DP = DT
               DT = FF - FO
               IF (DT.GT.DN) DN = DT
   80          CONTINUE
               FO = FN
  100       CONTINUE
            IF (BINOM) THEN
               IF2 = 0
               CALL G01BJF(NBIN,PAR(2),KPREV,FF,PGTK,PEQK,IF2)
               DT = FO - FF
               IF (DT.GT.DP) DP = DT
               DT = FF - FO
               IF (DT.GT.DN) DN = DT
            ELSE IF (POISON) THEN
               IF2 = 0
               CALL G01BKF(PAR(1),KPREV,FF,PGTK,PEQK,IF2)
               DT = FO - FF
               IF (DT.GT.DP) DP = DT
               DT = FF - FO
               IF (DT.GT.DN) DN = DT
            END IF
            DP = MAX(0.0D0,DP)
            DN = MAX(0.0D0,DN)
            IF (NTYPE.EQ.1) THEN
               D = MAX(DP,DN)
               P = 2.0D0*G08CBZ(N,D)
               IF (P.GT.1.0D0) P = 1.0D0
            ELSE IF (NTYPE.EQ.2) THEN
               D = DP
               P = G08CBZ(N,D)
            ELSE IF (NTYPE.EQ.3) THEN
               D = DN
               P = G08CBZ(N,D)
            END IF
            Z = SQRT(EN)*D
         END IF
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, N.lt.3 : N = ',I16)
99998 FORMAT (1X,'** On entry, DIST is not valid : DIST = ',A)
99997 FORMAT (1X,'** On entry, NTYPE is not equal to 1, 2 or 3: NTYPE ',
     *       '= ',I16)
99996 FORMAT (1X,'** On entry, ESTIMA is not valid : ESTIMA = ',A1)
99995 FORMAT (1X,'** On entry, the lower end of the range of the unifo',
     *       'rm distribution is not',/5X,'less than the upper end of ',
     *       'the range: a = PAR(1) = ',D13.5,/5X,'and  b = PAR(2) = ',
     *       D13.5,' .')
99994 FORMAT (1X,'** On entry, the variance for the normal dist. is le',
     *       'ss than zero',/5X,'PAR(2) = ',D13.5)
99993 FORMAT (1X,'** On entry, one of the gamma parameters is less or ',
     *       'equal to zero ',/5X,'PAR(1) = ',D13.5,' and PAR(2) = ',
     *       D13.5)
99992 FORMAT (1X,'** On entry, one of the beta parameters are not valid'
     *       ,/5X,'PAR(1) = ',D13.5,' and PAR(2) = ',D13.5)
99991 FORMAT (1X,'** On entry, n = PAR(1) for the binomial distributio',
     *       'n is not valid.',/5X,'Note that n must always be supplie',
     *       'd : PAR(1) = ',D13.5)
99990 FORMAT (1X,'** On entry, p = PAR(2) for the binomial distributio',
     *       'n is either less',/5X,'than or equal to zero or greater ',
     *       'than or equal to one',/5X,': PAR(2) = ',D13.5)
99989 FORMAT (1X,'** On entry, PAR(1) for the exponential distribution',
     *       ' is less than or ',/5X,'equal to zero : PAR(1) = ',D13.5)
99988 FORMAT (1X,'** On entry, PAR(1) for the Poisson distribution is ',
     *       'either less than or ',/5X,'equal to zero or greater than',
     *       ' 1000000 : PAR(1) = ',D13.5)
99987 FORMAT (1X,'** At least one of the sample observations is outsid',
     *       'e the range of the',/5X,'uniform distribution, i.e. less',
     *       ' than PAR(1) = ',D13.5,' or greater than',/5X,'PAR(2) = ',
     *       D13.5)
99986 FORMAT (1X,'** At least one of the sample observations is less t',
     *       'han zero.',/5X,'For testing against the gamma dist all o',
     *       'bservations should be non-negative.')
99985 FORMAT (1X,'** At least one of the sample observations is less t',
     *       'han zero or greater than',/5X,'one. For testing against ',
     *       'the beta distribution all observations should be',/5X,
     *       'in the (0,1) interval.')
99984 FORMAT (1X,'** At least one of the sample observations is less t',
     *       'han zero or greater than ',/5X,'n = ',D13.5,'. Remember ',
     *       'that for the binomial dist. you MUST provide a ',/5X,'va',
     *       'lue for n in the parameter PAR(1) and all observations m',
     *       'ust ',/5X,'be less than or equal to n.')
99983 FORMAT (1X,'** At least one of the sample observations is less t',
     *       'han zero. For testing',/5X,'against the exponential dist',
     *       '. all observations should be non-negative.')
99982 FORMAT (1X,'** At least one of the sample observations is less t',
     *       'han zero. For testing',/5X,'against the Poisson dist. al',
     *       'l observations should be non-negative.')
99981 FORMAT (1X,'** The whole sample is constant which implies the va',
     *       'riance is zero.',/5X,'This is only possible when using t',
     *       'he uniform, normal, gamma or beta dist.')
99980 FORMAT (1X,'** The variance = n*p*(1-p) for the binomial dist. e',
     *       'xceeds 1000000.',/5X,'n = PAR(1) = ',D13.5,' and p = PAR',
     *       '(2) = ',D13.5,' and VAR = ',D13.5)
99979 FORMAT (1X,'** In the computation of the incomplete gamma functi',
     *       'on by S14BAF the',/5X,'convergence of the Taylor series ',
     *       'or Legendre continued fraction',/5X,'fails within 600 it',
     *       'erations.')
99978 FORMAT (1X,'** All observations equal zero. When estimating para',
     *       'meters for the binomial',/'    dist., at least one obser',
     *       'vation must be greater than zero')
99977 FORMAT (1X,'** All observations equal zero. When estimating para',
     *       'meters for the exponential',/'    dist., at least one ob',
     *       'servation must be greater than zero')
99976 FORMAT (1X,'** All observations equal zero. When estimating para',
     *       'meters for the Poisson',/'    dist., at least one observ',
     *       'ation must be greater than zero')
      END
