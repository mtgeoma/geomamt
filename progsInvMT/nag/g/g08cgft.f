      SUBROUTINE G08CGF(K2,IFREQ,CINT,DIST,PAR,IPARAM,PROB,CHISQ,P,NDF,
     *                  EVAL,CHISQI,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08CGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHISQ, P
      INTEGER           IFAIL, IPARAM, K2, NDF
      CHARACTER         DIST
C     .. Array Arguments ..
      DOUBLE PRECISION  CHISQI(K2), CINT(K2-1), EVAL(K2), PAR(2),
     *                  PROB(K2)
      INTEGER           IFREQ(K2)
C     .. Local Scalars ..
      DOUBLE PRECISION  DF, STORP, STORP1, SUM, ULEN, XN
      INTEGER           I, IERR, IF2, IFAULT, N, NREC, NSMALL
      LOGICAL           CHI, EXPON, GAMMA, NORM, UNIF, USERDF
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF, G08CGZ, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, G08CGZ, X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
C
      CHISQ = 0.0D0
      IERR = 0
      NREC = 1
      N = 0
      NORM = .FALSE.
      UNIF = .FALSE.
      EXPON = .FALSE.
      CHI = .FALSE.
      GAMMA = .FALSE.
      USERDF = .FALSE.
C
C     Check input parameters.
C
      IF (K2.LT.2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) K2
      ELSE IF (DIST.NE.'N' .AND. DIST.NE.'n' .AND. DIST.NE.'U' .AND.
     *         DIST.NE.'u' .AND. DIST.NE.'E' .AND. DIST.NE.'e' .AND.
     *         DIST.NE.'C' .AND. DIST.NE.'c' .AND. DIST.NE.'G' .AND.
     *         DIST.NE.'g' .AND. DIST.NE.'A' .AND. DIST.NE.'a') THEN
         IERR = 2
         WRITE (P01REC,FMT=99998) DIST
      ELSE IF (IPARAM.LT.0 .OR. IPARAM.GE.K2-1) THEN
         NREC = 1
         IERR = 3
         WRITE (P01REC,FMT=99997) IPARAM
C
C        Check for elements of IFREQ.lt.0, and that CINT is in ascending
C          order, (and calculate N).
C
      ELSE
C
         IF (DIST.EQ.'N' .OR. DIST.EQ.'n') NORM = .TRUE.
         IF (DIST.EQ.'U' .OR. DIST.EQ.'u') UNIF = .TRUE.
         IF (DIST.EQ.'E' .OR. DIST.EQ.'e') EXPON = .TRUE.
         IF (DIST.EQ.'C' .OR. DIST.EQ.'c') CHI = .TRUE.
         IF (DIST.EQ.'G' .OR. DIST.EQ.'g') GAMMA = .TRUE.
         IF (DIST.EQ.'A' .OR. DIST.EQ.'a') USERDF = .TRUE.
C
         IF (IFREQ(1).LT.0) THEN
            IERR = 4
            WRITE (P01REC,FMT=99996)
            GO TO 140
         END IF
         N = N + IFREQ(1)
         DO 20 I = 2, K2 - 1
            N = N + IFREQ(I)
            IF (IFREQ(I).LT.0) THEN
               IERR = 4
               WRITE (P01REC,FMT=99996)
               GO TO 140
            END IF
            IF (CINT(I).LE.CINT(I-1)) THEN
               IERR = 5
               WRITE (P01REC,FMT=99995)
               GO TO 140
            END IF
   20    CONTINUE
         IF (IFREQ(K2).LT.0) THEN
            IERR = 4
            WRITE (P01REC,FMT=99996)
            GO TO 140
         END IF
         N = N + IFREQ(K2)
         IF ((EXPON .OR. CHI .OR. GAMMA) .AND. CINT(1).LT.0.0D0) THEN
            IERR = 6
            WRITE (P01REC,FMT=99994)
            GO TO 140
         END IF
C
C        Check the parameters of the CDF's
C
         NREC = 2
         IF (NORM .AND. PAR(2).LE.0.0D0) THEN
            IERR = 7
            WRITE (P01REC,FMT=99993) PAR(2)
         ELSE IF (UNIF) THEN
            IF (PAR(1).GE.PAR(2) .OR. PAR(1).GT.CINT(1) .OR. PAR(2)
     *          .LT.CINT(K2-1)) THEN
               IERR = 7
               WRITE (P01REC,FMT=99992) PAR(1), PAR(2)
            END IF
         ELSE IF (EXPON .AND. PAR(1).LE.0.0D0) THEN
            IERR = 7
            WRITE (P01REC,FMT=99991) PAR(1)
         ELSE IF (CHI .AND. PAR(1).LE.0.0D0) THEN
            IERR = 7
            WRITE (P01REC,FMT=99990) PAR(1)
         ELSE IF (GAMMA .AND. (PAR(1).LE.0.0D0 .OR. PAR(2).LE.0.0D0))
     *            THEN
            IERR = 7
            WRITE (P01REC,FMT=99989) PAR(1), PAR(2)
         ELSE IF (USERDF) THEN
            SUM = 0.0D0
            DO 40 I = 1, K2
               IF (PROB(I).LE.0.0D0) THEN
                  IERR = 8
                  WRITE (P01REC,FMT=99988) I, PROB(I)
                  GO TO 140
               END IF
               SUM = SUM + PROB(I)
   40       CONTINUE
            IF (ABS(SUM-1.0D0).GT.X02AJF()) THEN
               IERR = 8
               WRITE (P01REC,FMT=99987) SUM
               GO TO 140
            END IF
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         NREC = 0
         XN = DBLE(N)
C
C        Determine which distribution is to be used and calculate
C          probabilities, expected frequencies.
C
         IF (USERDF) THEN
C
C           User supplied probabilities.
C
            DO 60 I = 1, K2
               EVAL(I) = PROB(I)*XN
   60       CONTINUE
         ELSE IF (UNIF) THEN
C
C           Uniform distribution.
C
            ULEN = PAR(2) - PAR(1)
            EVAL(1) = XN*(CINT(1)-PAR(1))/ULEN
            DO 80 I = 2, K2 - 1
               EVAL(I) = XN*(CINT(I)-CINT(I-1))/ULEN
   80       CONTINUE
            EVAL(K2) = XN*(PAR(2)-CINT(K2-1))/ULEN
         ELSE
C
C           One of the other standard distributions.
C
            IFAULT = 0
            STORP1 = 0.0D0
            DO 100 I = 1, K2 - 1
               STORP = G08CGZ(DIST,CINT(I),PAR,IFAULT)
               IF (IFAULT.NE.0) THEN
                  NREC = 2
                  IERR = 11
                  WRITE (P01REC,FMT=99986)
               END IF
               EVAL(I) = (STORP-STORP1)*XN
               STORP1 = STORP
  100       CONTINUE
            STORP = 1.0D0 - STORP
            EVAL(K2) = STORP*XN
         END IF
C
C        Calculate contributions to the statistic and hence
C        the statistic.
C
         NSMALL = 0
         DO 120 I = 1, K2
            IF (EVAL(I).LT.1.0D0) NSMALL = NSMALL + 1
            IF (EVAL(I).NE.0.0D0) THEN
               CHISQI(I) = (DBLE(IFREQ(I))-EVAL(I))*(DBLE(IFREQ(I))
     *                     -EVAL(I))/EVAL(I)
            ELSE IF (EVAL(I).EQ.0.0D0 .AND. IFREQ(I).EQ.0) THEN
               CHISQI(I) = 0.0D0
            ELSE
               IERR = 9
               NREC = 1
               WRITE (P01REC,FMT=99985)
               GO TO 140
            END IF
            CHISQ = CHISQ + CHISQI(I)
  120    CONTINUE
         IF (NSMALL.NE.0) THEN
            IERR = 10
            NREC = 1
            WRITE (P01REC,FMT=99984) NSMALL
         END IF
C
C        Calculate the degrees of freedom associated with the test.
C        No. intervals - 1 - No. estimated parameters.
C
         NDF = K2 - 1 - IPARAM
         DF = DBLE(NDF)
C
C        Determine the probability associated with the statistic.
C
         IF2 = 1
         P = G01ECF('UPPER',CHISQ,DF,IF2)
      END IF
  140 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, NCLASS.lt.2 : NCLASS = ',I16)
99998 FORMAT (1X,'** On entry, DIST is not valid : DIST = ',A1)
99997 FORMAT (1X,'** On entry, NPEST.lt.0 or NPEST.ge.NCLASS-1 : NPEST',
     *       ' = ',I16)
99996 FORMAT (1X,'** On entry, at least one class has frequency less t',
     *       'han zero.')
99995 FORMAT (1X,'** On entry, the elements of CINT are not in ascendi',
     *       'ng order.')
99994 FORMAT (1X,'** On entry, the intervals contained in CINT are inv',
     *       'alid.')
99993 FORMAT (1X,'** On entry, the variance of the normal distribution',
     *       ',PAR(2) is invalid:',/3X,'PAR(2) = ',1P,D13.5)
99992 FORMAT (1X,'** On entry, the parameters of the uniform  distribu',
     *       'tion are invalid:',/3X,'PAR(1) = ',1P,D13.5,' PAR(2) = ',
     *       D13.5)
99991 FORMAT (1X,'** On entry, the parameter of the exponential distri',
     *       'bution is invalid:',/4X,'PAR(1) = ',1P,D13.5)
99990 FORMAT (1X,'** On entry, the parameter of the Chi-squared distri',
     *       'bution is invalid:',/4X,'PAR(1) = ',1P,D13.5)
99989 FORMAT (1X,'** On entry, the parameter(s) of the gamma distribut',
     *       'ion are invalid: ',/4X,'PAR(1) = ',1P,D13.5,' PAR(2) = ',
     *       D13.5)
99988 FORMAT (1X,'** On entry, with DIST.eq.''A'' or ''a'' at least on',
     *       'e elements of PROB.le.0.0',/4X,'For the ',I16,' th class',
     *       ' PROB(i) = ',D13.5)
99987 FORMAT (1X,'** On entry with DIST.eq.''A'' or ''a'' the sum of t',
     *       'he elements of PROB',/4X,'does .ne. 1.0.  SUM of PROB(i)',
     *       '  = ',D13.5)
99986 FORMAT (1X,'** The solution has failed to converge whilst comput',
     *       'ing expected values for the',/4X,'gamma orchi-squared di',
     *       'st. The result returned should be an adequate approx.')
99985 FORMAT (1X,'** An expected frequency equals zero, when the obser',
     *       'ved frequency was not.')
99984 FORMAT (1X,'** ',I16,' classes have expected frequency less than',
     *       ' one.')
      END
