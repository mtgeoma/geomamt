      SUBROUTINE G07EBF(METHOD,N,X,M,Y,CLEVEL,THETA,THETAL,THETAU,ESTCL,
     *                  ULOWER,UUPPER,WRK,IWRK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07EBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CLEVEL, ESTCL, THETA, THETAL, THETAU, ULOWER,
     *                  UUPPER
      INTEGER           IFAIL, M, N
      CHARACTER         METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(3*(N+M)), X(N), Y(M)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALNEW, ALPHA, AUNEW, D, DF, DH, DL, EPS, PROBL,
     *                  PROBU, SD, SLOPE, T, THETA2, TOL, TOL1, TOL2,
     *                  TPROB, U, UNOR, VARXB, VARYB, WIDTH, XBAR, XMP,
     *                  YBAR
      INTEGER           I, IERROR, IF2, IFAIL2, ILOWER, IND, IND1, IND2,
     *                  IUPPER, LWRK, MMINUS, NMINUS, NREC
      LOGICAL           EXACT, LIMITS, LOWER, TIES
C     .. Local Arrays ..
      DOUBLE PRECISION  WRK2(201)
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  G01FBF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01FBF, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06EFF, G07EBT, G07EBW, G07EBY, G07EBZ, G08AHF,
     *                  G08AJF, M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, MOD, SQRT
C     .. Executable Statements ..
      THETA = 0.0D0
      THETAL = 0.0D0
      THETAU = 0.0D0
C
C     Check for invalid character, insufficient data, or invalid
C     percent confidence.
C
      IERROR = 0
      NREC = 0
      IF (METHOD.NE.'E' .AND. METHOD.NE.'e' .AND. METHOD.NE.'A' .AND.
     *    METHOD.NE.'a') THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) METHOD
      ELSE IF (N.LT.1 .OR. M.LT.1) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99998) N, M
      ELSE IF (CLEVEL.GE.1.0D0 .OR. CLEVEL.LE.0.0D0) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99997) CLEVEL
      ELSE IF (N.EQ.1 .AND. M.EQ.1) THEN
         THETA = Y(1) - X(1)
         THETAL = THETA
         THETAU = THETA
         ULOWER = 0.0D0
         UUPPER = 1.0D0
         ESTCL = 0.0D0
      ELSE IF (N.EQ.1 .AND. M.EQ.2) THEN
         THETA = (Y(1)+Y(2))/2.0D0 - X(1)
         THETAL = MIN(Y(1),Y(2)) - X(1)
         THETAU = MAX(Y(1),Y(2)) - X(1)
         ULOWER = 0.0D0
         UUPPER = 2.0D0
         ESTCL = 1.0D0/3.0D0
      ELSE IF (N.EQ.2 .AND. M.EQ.1) THEN
         THETA = Y(1) - (X(1)+X(2))/2.0D0
         THETAL = Y(1) - MAX(X(1),X(2))
         THETAU = Y(1) - MIN(X(1),X(2))
         ULOWER = 0.0D0
         UUPPER = 2.0D0
         ESTCL = 1.0D0/3.0D0
      ELSE
         CALL F06EFF(N,X,1,WRK,1)
         IF2 = 0
         CALL M01CAF(WRK,1,N,'A',IF2)
         CALL F06EFF(M,Y,1,WRK(N+1),1)
         CALL M01CAF(WRK,N+1,N+M,'A',IF2)
         IF (WRK(1).EQ.WRK(N) .AND. WRK(N+1).EQ.WRK(N+M)) THEN
C
C           All values in the first sample identical and all values in
C           the second sample identical.
C
            THETA = WRK(N+1) - WRK(1)
            THETAL = THETA
            THETAU = THETA
            ULOWER = 0.0D0
            UUPPER = DBLE(N*M)
            ESTCL = 0.0D0
            IERROR = 2
            NREC = 3
            WRITE (P01REC,FMT=99996)
         ELSE
C
C           Preliminary estimates using trimmed means and variances
C
            CALL G07EBT(WRK,N,XBAR,VARXB)
            CALL G07EBT(WRK(N+1),M,YBAR,VARYB)
            D = YBAR - XBAR
            SD = SQRT(VARXB+VARYB)
            SD = MAX(SD,1.0D-20)
            DF = DBLE(M+N-2)
C
C           ALPHA is one tailed probability.
C
            ALPHA = 0.5D0 - CLEVEL/2.0D0
            IF2 = 0
            T = G01FBF('Upper',ALPHA,DF,IF2)
            WIDTH = T*SD
            IF (WIDTH/D.GT.SQRT(X02AJF())) THEN
               DL = D - WIDTH
               DH = D + WIDTH
            ELSE
               NMINUS = MAX(1,N/2-1)
               MMINUS = MAX(1,M/2-1)
               DH = WRK(N+M/2+1) - WRK(NMINUS)
               DL = WRK(N+MMINUS) - WRK(N/2+1)
               D = (DH+DL)/2.0D0
            END IF
C
C           Define tolerance TOL1 based on the magnitude of the data
C           Define tolerance TOL2 based on width of confidence interval
C           (based on preliminary estimates). The point estimate should
C           be accurate to 5 significant digits in the width (approx.).
C           Maximum accuracy set to 1.E-16
C
C
            XMP = X02AJF()
            TOL1 = 2.0D0*XMP*MAX(D,DL,DH)
            TOL2 = 2.0D0*SD*0.25D-05
            TOL = MAX(TOL1,TOL2,1.0D-16)
            EPS = 2.0D0*XMP
C
C           Find the values of the Mann-Whitney U test statistic
C           satisfying required tail probabilities.
C
            EXACT = (N+M) .LE. 40 .OR. MAX(N,M) .LE. 30
            IF (EXACT) THEN
               LWRK = N*M/2 + 1
            ELSE
               LWRK = 1
            END IF
            LOWER = .TRUE.
            CALL G07EBZ(LOWER,ALPHA,N,M,ILOWER,WRK2,LWRK)
            ULOWER = DBLE(ILOWER)
            LOWER = .FALSE.
            CALL G07EBZ(LOWER,ALPHA,N,M,IUPPER,WRK2,LWRK)
            UUPPER = DBLE(IUPPER)
C
C           METHOD = 'E' - use exact algorithm
C           METHOD = 'A' - use root finding approximation
C
            IF (METHOD.EQ.'E' .OR. METHOD.EQ.'e') THEN
C
C              Find the lower confidence limit, THETAL.
C
               IND = ILOWER + 1
               CALL G07EBW(N,WRK,M,WRK(N+1),IND,IND,THETAL,IWRK)
C
C              Find the upper confidence limit, THETAU.
C
               IND = IUPPER
               CALL G07EBW(N,WRK,M,WRK(N+1),IND,IND,THETAU,IWRK)
C
C              Now find the estimate of location THETA.
C
               IF (MOD((N*M),2).EQ.0) THEN
                  IND1 = N*M/2
                  IND2 = IND1 + 1
               ELSE
                  IND1 = (N*M+1)/2
                  IND2 = IND1
               END IF
               CALL G07EBW(N,WRK,M,WRK(N+1),IND1,IND2,THETA,IWRK)
            ELSE
C
C              ULOWER  and UUPPER  give the ordered differences which
C              form the confidence interval. An estimate of the slope
C              of the linear approximation to G07EBY.
C
               SLOPE = (UUPPER-ULOWER)/(DH-DL)
C
C              Find the lower confidence limit, THETAL.
C
               IND = ILOWER
               LIMITS = .TRUE.
               LOWER = .TRUE.
               CALL G07EBY(N,WRK,M,WRK(N+1),IND,DL,LIMITS,LOWER,SLOPE,
     *                     THETAL,EPS,TOL,IERROR)
               IF (IERROR.EQ.3) THEN
                  NREC = 2
                  WRITE (P01REC,FMT=99995)
               END IF
C
C              Find the upper confidence limit, THETAU.
C
               IND = IUPPER
               LOWER = .FALSE.
               CALL G07EBY(N,WRK,M,WRK(N+1),IND,DH,LIMITS,LOWER,SLOPE,
     *                     THETAU,EPS,TOL,IERROR)
               IF (IERROR.EQ.3) THEN
                  NREC = 2
                  WRITE (P01REC,FMT=99994)
               END IF
C
C              A new estimate of slope based on the confidence interval
C              (THETAL,THETAU), unless the length of the interval is
C              smaller than EPS.
C
               IF (THETAU.GT.THETAL+EPS) SLOPE = ((DH-DL)
     *             /(THETAU-THETAL))*SLOPE
C
C              The mid-point of the confidence interval will be the
C              initial estimate of THETA.
C
               D = (THETAL+THETAU)/2.0D0
               LIMITS = .FALSE.
               LOWER = .FALSE.
               IND = (N*M+1)/2
               CALL G07EBY(N,WRK,M,WRK(N+1),IND,D,LIMITS,LOWER,SLOPE,
     *                     THETA,EPS,TOL,IERROR)
               IF (IERROR.EQ.3) THEN
                  NREC = 2
                  WRITE (P01REC,FMT=99993)
               END IF
C
               IF (MOD(N*M,2).EQ.0) THEN
                  LOWER = .TRUE.
                  CALL G07EBY(N,WRK,M,WRK(N+1),IND,D,LIMITS,LOWER,SLOPE,
     *                        THETA2,EPS,TOL,IERROR)
                  IF (IERROR.EQ.3) THEN
                     NREC = 2
                     WRITE (P01REC,FMT=99993)
                  END IF
                  THETA = (THETA+THETA2)/2.0D0
               END IF
C
            END IF
C
C           Now estimate the percentage confidence for the interval
C           found.
C
            ALNEW = THETAL - 10.0D0*XMP
            DO 20 I = 1, N
               WRK(I) = WRK(I) + ALNEW
   20       CONTINUE
            CALL G08AHF(N,WRK,M,WRK(N+1),'Lower',U,UNOR,PROBL,TIES,
     *                  WRK(N+M+1),WRK(2*(N+M)+1),IFAIL2)
            IF (EXACT) CALL G08AJF(N,M,'Lower',U,PROBL,WRK(N+M+1),LWRK,
     *                             IFAIL2)
C
            AUNEW = THETAU + 10.0D0*XMP
            DO 40 I = 1, N
               WRK(I) = X(I) + AUNEW
   40       CONTINUE
            CALL G08AHF(N,WRK,M,WRK(N+1),'Upper',U,UNOR,PROBU,TIES,
     *                  WRK(N+M+1),WRK(2*(N+M)+1),IFAIL2)
            IF (EXACT) CALL G08AJF(N,M,'Upper',U,PROBU,WRK(N+M+1),LWRK,
     *                             IFAIL2)
            TPROB = PROBL + PROBU
            ESTCL = 1.0D0 - TPROB
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, METHOD is not valid: METHOD = ',A1)
99998 FORMAT (1X,'** On entry, N.lt.1 or M.lt.1 : N = ',I16,', M = ',
     *       I16)
99997 FORMAT (1X,'** On entry, CLEVEL is out of range : CLEVEL = ',
     *       D13.5)
99996 FORMAT (1X,'** Not enough information to compute an interval est',
     *       'imate since each',/4X,'sample has identical values. The ',
     *       'common difference is returned in',/4X,'THETA, THETAL and',
     *       ' THETAU.')
99995 FORMAT (1X,'** Warning. The iterative procedure to find an estim',
     *       'ate of the lower ',/4X,'confidence limit has not converg',
     *       'ed in 100 iterations.')
99994 FORMAT (1X,'** Warning. The iterative procedure to find an estim',
     *       'ate of the upper ',/4X,'confidence limit has not converg',
     *       'ed in 100 iterations.')
99993 FORMAT (1X,'** Warning. The iterative procedure to find an estim',
     *       'ate of Theta ',/4X,'has not converged in 100 iterations.')
      END
