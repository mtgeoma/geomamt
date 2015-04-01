      SUBROUTINE G13DJF(K,N,Z,IK,TR,ID,DELTA,IP,IQ,MEAN,PAR,LPAR,QQ,V,
     *                  LMAX,PREDZ,SEFZ,REF,LREF,WORK,LWORK,IWORK,
     *                  LIWORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 15B REVISED. IER-957 (NOV 1991).
C     MARK 17 REVISED. IER-1687 (JUN 1995).

C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, IP, IQ, K, LIWORK, LMAX, LPAR, LREF,
     *                  LWORK, N
      CHARACTER         MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  DELTA(IK,*), PAR(LPAR), PREDZ(IK,LMAX),
     *                  QQ(IK,K), REF(LREF), SEFZ(IK,LMAX), V(IK,*),
     *                  WORK(LWORK), Z(IK,N)
      INTEGER           ID(K), IWORK(LIWORK)
      CHARACTER         TR(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGEST, LIMIT, SUM2
      INTEGER           I, IDMAX, IERR, IERROR, IFAILA, IR, J, JF, JREF,
     *                  KP, KQ, LR, LR1, LW1, LW2, LW3, LW4, ND, NOBS,
     *                  NP, NPAR, NREC
      LOGICAL           STAT
      CHARACTER         AA
C     .. Local Arrays ..
      CHARACTER*80      P01REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03ABF, DCOPY, G05HDZ, G13DJZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     This subroutine calculates linear minimum mean square error
C     forecasts of future series values using the parameter estimates
C     and the residuals returned by G13DCF. The standard deviations
C     of the forecast errors are also returned and a reference vector
C     (REF) set up for G13DKF.
C
C     First test for errors in the input arguments
C
      IERROR = 1
      NREC = 1
      NPAR = (IP+IQ)*K*K
      IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') NPAR = NPAR + K
      NP = NPAR + K*(K+1)/2
      NOBS = N*K
      IF (K.LT.1) THEN
         WRITE (P01REC,FMT=99999) K
      ELSE IF (N.LT.3) THEN
         WRITE (P01REC,FMT=99998) N
      ELSE IF (IK.LT.K) THEN
         WRITE (P01REC,FMT=99997) IK, K
      ELSE IF (IP.LT.0) THEN
         WRITE (P01REC,FMT=99996) IP
      ELSE IF (IQ.LT.0) THEN
         WRITE (P01REC,FMT=99995) IQ
      ELSE IF (MEAN.NE.'M' .AND. MEAN.NE.'m' .AND. MEAN.NE.'Z' .AND.
     *         MEAN.NE.'z') THEN
         WRITE (P01REC,FMT=99994)
      ELSE IF (LPAR.LT.NPAR) THEN
         NREC = 2
         WRITE (P01REC,FMT=99993) LPAR, NPAR
      ELSE IF (NOBS.LE.NP) THEN
         NREC = 3
         WRITE (P01REC,FMT=99992) NOBS, NP
      ELSE IF (LMAX.LT.1) THEN
         WRITE (P01REC,FMT=99991) LMAX
      ELSE
         IERROR = 0
         NREC = 0
         IR = MAX(IP,IQ)
         IDMAX = 0
         JF = (LMAX-1)*K*K + 2*K*LMAX
C
C        test whether LREF is big enough
C
         JREF = JF + K
         IF (LREF.LT.JREF) THEN
            IERROR = 1
            NREC = 2
            WRITE (P01REC,FMT=99987) LREF, JREF
            GO TO 80
         END IF
C
         DO 20 I = 1, K
C
C           Check that all elements of ID are .ge. 0 and less than
C           N - max(IP,IQ) and find the maximum element of ID.
C
            IF (ID(I).LT.0) THEN
               IERROR = 1
               NREC = 1
               WRITE (P01REC,FMT=99990) I
               GO TO 80
            ELSE IF (IR+ID(I).GE.N) THEN
               IERROR = 1
               NREC = 1
               WRITE (P01REC,FMT=99989) I
               GO TO 80
            ELSE
               IDMAX = MAX(ID(I),IDMAX)
            END IF
C
C           Test for valid TR array
C
            AA = TR(I)
            IF (AA.EQ.'N' .OR. AA.EQ.'n') THEN
               REF(JF+I) = 100.0D0
            ELSE IF (AA.EQ.'L' .OR. AA.EQ.'l') THEN
               REF(JF+I) = 200.0D0
            ELSE IF (AA.EQ.'S' .OR. AA.EQ.'s') THEN
               REF(JF+I) = 300.0D0
            ELSE
               IERROR = 2
               NREC = 1
               WRITE (P01REC,FMT=99988) I
               GO TO 80
            END IF
   20    CONTINUE
C
C        test whether real and integer workspace arrays are big enough
C
         I = IR
         I = ((IR*K)**2) + 2*K*IR
         J = ((IP+IDMAX+2)*K*K) + ((N+LMAX)*K)
         J = MAX(I,J)
C
         IF (LWORK.LT.J) THEN
            IERROR = 1
            NREC = 2
            WRITE (P01REC,FMT=99986) LWORK, J
            GO TO 80
         END IF
C
         J = K*IR
C
         IF (LIWORK.LT.J) THEN
            IERROR = 1
            NREC = 2
            WRITE (P01REC,FMT=99985) LIWORK, J
            GO TO 80
         END IF
C
         ND = N - IDMAX
C
C        test whether QQ is positive-definite
C
C        first set the upper triangle of QQ to the lower triangle of QQ
C
         DO 40 I = 1, K - 1
            CALL DCOPY(K-I,QQ(I+1,I),1,QQ(I,I+1),IK)
   40    CONTINUE
C
         IFAILA = 1
         CALL F03ABF(QQ,IK,K,SUM2,WORK,IFAILA)
C
C        reset the strict lower triangle of QQ (which F03ABF has
C        corrupted) to its former value
C
         DO 60 I = 1, K - 1
            CALL DCOPY(K-I,QQ(I,I+1),IK,QQ(I+1,I),1)
   60    CONTINUE
C
         IF (IFAILA.GT.0) THEN
            IERROR = 4
            NREC = 1
            WRITE (P01REC,FMT=99984)
            GO TO 80
         END IF
C
C        test for stationarity of the AR operator
C
         IERR = 0
         IF (IP.GT.0) THEN
            KP = K*IP
            LW1 = 1
            LW2 = LW1 + KP*KP
            LW3 = LW2 + KP
            LIMIT = 1.0D0
            CALL G05HDZ(IP,K,PAR,WORK(LW1),WORK(LW2),WORK(LW3),IWORK,KP,
     *                  LIMIT,BIGEST,STAT,IERR)
            IF (IERR.EQ.1) THEN
               IERROR = 5
               NREC = 2
               WRITE (P01REC,FMT=99983)
               GO TO 80
            ELSE IF ( .NOT. STAT) THEN
               IERROR = 4
               NREC = 1
               WRITE (P01REC,FMT=99982)
               GO TO 80
            END IF
         END IF
C
C        test for invertibility of the MA operator
C
         IF (IQ.GT.0) THEN
            LW4 = IP*K*K + 1
            KQ = K*IQ
            LW1 = 1
            LW2 = LW1 + KQ*KQ
            LW3 = LW2 + KQ
            LIMIT = 1.0D0
            CALL G05HDZ(IQ,K,PAR(LW4),WORK(LW1),WORK(LW2),WORK(LW3),
     *                  IWORK,KQ,LIMIT,BIGEST,STAT,IERR)
            IF (IERR.EQ.1) THEN
               IERROR = 5
               NREC = 2
               WRITE (P01REC,FMT=99983)
               GO TO 80
            ELSE IF ( .NOT. STAT) THEN
               IERROR = 4
               NREC = 1
               WRITE (P01REC,FMT=99981)
               GO TO 80
            END IF
         END IF
C
         LW1 = K*(N+LMAX) + 1
         LW2 = LW1 + K*K
         LW3 = LW2
         IF ((IP+IDMAX).NE.0) LW3 = LW3 + K*K
         LR1 = (LMAX-1)*K*K + 1
         LR = 2*LMAX*K
C
         CALL G13DJZ(K,N,IP,IQ,MEAN,PAR,NPAR,QQ,IK,Z,TR,ID,DELTA,IDMAX,
     *               V,LMAX,ND,PREDZ,SEFZ,WORK,WORK(LW1),WORK(LW2),
     *               WORK(LW3),REF,REF(LR1),LR,IERR)
C
         IF (IERR.EQ.1) THEN
            IERROR = 3
            NREC = 1
            WRITE (P01REC,FMT=99980)
         ELSE IF (IERR.EQ.2) THEN
            IERROR = 6
            NREC = 1
            WRITE (P01REC,FMT=99979)
         ELSE IF (IERR.EQ.3) THEN
            IERROR = 7
            NREC = 1
            WRITE (P01REC,FMT=99978)
         END IF
      END IF
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT ('  ** On entry, K.lt.1 : K = ',I16)
99998 FORMAT ('  ** On entry, N.lt.3 : N = ',I16)
99997 FORMAT ('  ** On entry, IK.lt.K : IK = ',I16,'  and K = ',I16)
99996 FORMAT ('  ** On entry, IP.lt.0 : IP = ',I16)
99995 FORMAT ('  ** On entry, IQ.lt.0 : IQ = ',I16)
99994 FORMAT ('  ** On entry, MEAN is an invalid character.')
99993 FORMAT ('  ** On entry, LPAR is too small : LPAR = ',I16,' but m',
     *       'ust be at',/'     least ',I16)
99992 FORMAT ('  ** On entry, the total number of observations is less',
     *       ' than the total number of',/'    parameters (including t',
     *       'he covariance matrix).',/'    No. of obs. = ',I16,'  and',
     *       ' no. of parameters = ',I16)
99991 FORMAT (' ** On entry, LMAX.lt.1 : LMAX = ',I16)
99990 FORMAT ('  ** On entry, element ',I3,' of ID is less than zero')
99989 FORMAT ('  ** On entry, element ',I3,' of ID is greater than or ',
     *       'equal to N-max(IP,IQ)')
99988 FORMAT ('  ** On entry, element ',I3,' of the array TR is an inv',
     *       'alid character')
99987 FORMAT ('  ** On entry, LREF is too small : LREF = ',I16,'  but ',
     *       'must be at',/'     least ',I16)
99986 FORMAT ('  ** On entry, LWORK is too small : LWORK = ',I16,'  bu',
     *       't must be at',/'     least ',I16)
99985 FORMAT ('  ** On entry, LIWORK is too small : LIWORK = ',I16,'  ',
     *       'but must be',/'     at least ',I16)
99984 FORMAT ('  ** On entry, the covariance matrix QQ is not positive',
     *       '-definite.')
99983 FORMAT ('  ** An excessive number of iterations were needed by F',
     *       '02AFF to evaluate the',/'    eigenvalues of the matrices',
     *       ' used to test for stationarity and invertibility.')
99982 FORMAT ('  ** On entry, the AR parameter matrices are outside th',
     *       'e stationarity region.')
99981 FORMAT ('  ** On entry, the MA parameter matrices are outside th',
     *       'e invertibility region.')
99980 FORMAT ('  ** On entry, one (or more) of the transformations req',
     *       'uested is invalid.')
99979 FORMAT ('  ** The covariance matrix may be nearly non positive-d',
     *       'efinite.')
99978 FORMAT ('  ** The forecasts will overflow if computed.')
      END
