      SUBROUTINE G07DCF(CHI,PSI,ISIGMA,N,X,BETA,THETA,SIGMA,MAXIT,TOL,
     *                  RS,NIT,WRK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 13B REVISED. IER-666 (AUG 1988).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07DCF')
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.00D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BETA, SIGMA, THETA, TOL
      INTEGER           IFAIL, ISIGMA, MAXIT, N, NIT
C     .. Array Arguments ..
      DOUBLE PRECISION  RS(N), WRK(N), X(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  CHI, PSI
      EXTERNAL          CHI, PSI
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST, D, DSIG, RSS, S2, SCHI, SIGMB, TOL1, XMD,
     *                  XME, XSD
      INTEGER           I, IERROR, IFAIL2, NREC
      LOGICAL           ALLZER
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G07DAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, DBLE, SQRT
C     .. Executable Statements ..
      NREC = 1
C
C     PARAMETER CHECK
C
      IERROR = 1
      IF (N.LT.2) THEN
         WRITE (REC,FMT=99999) N
      ELSE IF (TOL.LE.ZERO) THEN
         WRITE (REC,FMT=99998) TOL
      ELSE IF (MAXIT.LE.0) THEN
         WRITE (REC,FMT=99997) MAXIT
      ELSE IF (ISIGMA.NE.0 .AND. ISIGMA.NE.1) THEN
         WRITE (REC,FMT=99996) ISIGMA
      ELSE IF (BETA.LE.ZERO) THEN
         IERROR = 2
         WRITE (REC,FMT=99991) BETA
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
C
C        CHECK FOR ALL ELEMENTS OF X BEING EQUAL
C
         DO 20 I = 2, N
            IF (X(I).NE.X(I-1)) GO TO 40
   20    CONTINUE
         IERROR = 3
         WRITE (REC,FMT=99994)
         GO TO 240
C
C        CALCULATE STARTING VALUES
C
   40    IFAIL2 = 0
         IF (SIGMA.LE.ZERO) THEN
            CALL G07DAF(N,X,WRK,XME,XMD,XSD,IFAIL2)
C           THIS IFAIL CANNOT BE TRIPPED.
            THETA = XME
            SIGMA = XSD
            SIGMB = XSD
         ELSE
            SIGMB = SIGMA
         END IF
         CONST = BETA*DBLE(N-1)
C
C        STEP1 - SET NIT = 1
C
         NIT = 1
   60    CONTINUE
C
C        STEP2  -  COMPUTE RESIDUALS X-THETA
C
         DO 80 I = 1, N
            RS(I) = X(I) - THETA
   80    CONTINUE
C
C        STEP3  -  COMPUTE A NEW VALUE SIGMB FOR SIGMA
C
         IF (ISIGMA.NE.0) THEN
            S2 = ZERO
            DO 100 I = 1, N
               RSS = RS(I)/SIGMA
               SCHI = CHI(RSS)
               IF (SCHI.GE.ZERO) THEN
                  S2 = S2 + SCHI
               ELSE
                  IERROR = 7
                  WRITE (REC,FMT=99990) SCHI
                  GO TO 240
               END IF
  100       CONTINUE
            IF (S2.LE.ZERO) THEN
               GO TO 220
            ELSE
               SIGMB = SQRT(S2/CONST)*SIGMA
            END IF
         END IF
C
C        STEP4  -  WINSORIZE THE RESIDUALS
C
         DO 120 I = 1, N
            RSS = RS(I)/SIGMB
            RS(I) = PSI(RSS)*SIGMB
  120    CONTINUE
C
C        STEP5  -  COMPUTE THE INCREMENT
C
         D = ZERO
         DO 140 I = 1, N
            D = D + RS(I)
  140    CONTINUE
         D = D/DBLE(N)
C
C        STEP6  -  UPDATE THETA
C
         THETA = THETA + D
C
C        STEP7  -  STOP IF DESIRED PRECISION IS REACHED
C
         TOL1 = TOL*MAX(1.D0,SIGMB)
         IF (ABS(D).LT.TOL1) THEN
            DSIG = SIGMB - SIGMA
            IF (ABS(DSIG).LT.TOL1) GO TO 160
         END IF
C
C        ITERATIONS CONTINUE IF NIT. LT. MAXIT
C
         SIGMA = SIGMB
         IF (NIT.LT.MAXIT) THEN
            NIT = NIT + 1
            GO TO 60
         ELSE
            IERROR = 5
            NREC = 1
            WRITE (REC,FMT=99992) MAXIT
         END IF
C
C        ITERATIONS ARE TERMINATED
C
  160    CONTINUE
C
C        TEST THAT THE ITERATIONS WERE NOT TERMINATED
C        DUE TO ZERO RESIDUALS
C
         ALLZER = .TRUE.
         DO 180 I = 1, N
            IF (RS(I).NE.ZERO) ALLZER = .FALSE.
  180    CONTINUE
         IF (ALLZER) THEN
            IERROR = 6
            WRITE (REC,FMT=99993)
         END IF
         DO 200 I = 1, N
            RS(I) = X(I) - THETA
  200    CONTINUE
         SIGMA = SIGMB
         GO TO 240
  220    IERROR = 4
         WRITE (REC,FMT=99995) S2
      END IF
C
C     ASSIGN IFAIL VALUE
C
  240 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.le.1: N = ',I16)
99998 FORMAT (1X,'** On entry, TOL.le.0.0: TOL = ',1P,D13.5)
99997 FORMAT (1X,'** On entry, MAXIT.le.0: MAXIT = ',I16)
99996 FORMAT (1X,'** On entry, ISIGMA.ne.0 or 1: ISIGMA = ',I16)
99995 FORMAT (1X,'** Current estimate of SIGMA is zero or negative: SI',
     *       'GMA=',1P,D13.5)
99994 FORMAT (1X,'** All elements of X are equal.')
99993 FORMAT (1X,'** All winsorized residuals are zero.')
99992 FORMAT (1X,'** Number of iterations required exceeds MAXIT: MAXI',
     *       'T =',I12)
99991 FORMAT (1X,'** On entry, BETA.le.0.0: BETA = ',1P,D13.5)
99990 FORMAT (1X,'** The CHI function returned a negative value: CHI = '
     *       ,1P,D13.5)
      END
