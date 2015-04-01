      SUBROUTINE G07DBF(ISIGMA,N,X,IPSI,C,H1,H2,H3,DCHI,THETA,SIGMA,
     *                  MAXIT,TOL,RS,NIT,WRK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07DBF')
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.00D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, DCHI, H1, H2, H3, SIGMA, THETA, TOL
      INTEGER           IFAIL, IPSI, ISIGMA, MAXIT, N, NIT
C     .. Array Arguments ..
      DOUBLE PRECISION  RS(N), WRK(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, CONST, D, DSIG, RSS, S2, SIGMB, TOL1, XMD,
     *                  XME, XSD
      INTEGER           I, IERROR, IFAIL2, NREC
      LOGICAL           ALLZER
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G07DBX, G07DBY
      INTEGER           P01ABF
      EXTERNAL          G07DBX, G07DBY, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G07DAF, G07DBZ
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
      ELSE IF (IPSI.LT.0 .OR. IPSI.GT.4) THEN
         WRITE (REC,FMT=99995) IPSI
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
         WRITE (REC,FMT=99992)
         GO TO 240
   40    IERROR = 2
         IF (IPSI.EQ.1) THEN
            IF (C.LE.ZERO) THEN
               WRITE (REC,FMT=99994) C
               GO TO 240
            END IF
         END IF
         IF (IPSI.EQ.2) THEN
            IF ((H1.LT.ZERO .OR. (H1.EQ.H2 .AND. H1.EQ.H3 .AND. H1.EQ.
     *          ZERO) .OR. H1.GT.H2 .OR. H1.GT.H3 .OR. H2.GT.H3)) THEN
               NREC = 2
               WRITE (REC,FMT=99989) H1, H2, H3
               GO TO 240
            END IF
         END IF
         IF (IPSI.NE.0) THEN
            IF (DCHI.LE.ZERO) THEN
               WRITE (REC,FMT=99988) DCHI
               GO TO 240
            END IF
         END IF
         IERROR = 0
C
C        CALCULATE BETA
C
         IF (IERROR.EQ.0) THEN
            IF (IPSI.NE.0) THEN
               CALL G07DBZ(DCHI,BETA)
            ELSE
               BETA = 0.5D0
            END IF
C
C           CALCULATE STARTING VALUES
C
            IFAIL2 = 0
            IF (SIGMA.LE.ZERO) THEN
               CALL G07DAF(N,X,WRK,XME,XMD,XSD,IFAIL2)
C              THIS IFAIL CANNOT BE TRIPPED.
               THETA = XME
               SIGMA = XSD
               SIGMB = XSD
            ELSE
               SIGMB = SIGMA
            END IF
            CONST = BETA*DBLE(N-1)
C
C           STEP1 - SET NIT = 1
C
            NIT = 1
   60       CONTINUE
C
C           STEP2  -  COMPUTE RESIDUALS X-THETA
C
            DO 80 I = 1, N
               RS(I) = X(I) - THETA
   80       CONTINUE
C
C           STEP3  -  COMPUTE A NEW VALUE SIGMB FOR SIGMA
C
            IF (ISIGMA.NE.0) THEN
               S2 = ZERO
               DO 100 I = 1, N
                  RSS = RS(I)/SIGMA
                  S2 = S2 + G07DBY(RSS,IPSI,DCHI)
  100          CONTINUE
               IF (S2.LE.ZERO) THEN
                  GO TO 220
               ELSE
                  SIGMB = SQRT(S2/CONST)*SIGMA
               END IF
            END IF
C
C           STEP4  -  WINSORIZE THE RESIDUALS
C
            DO 120 I = 1, N
               RSS = RS(I)/SIGMB
               RS(I) = G07DBX(RSS,IPSI,C,H1,H2,H3)*SIGMB
  120       CONTINUE
C
C           STEP5  -  COMPUTE THE INCREMENT
C
            D = ZERO
            DO 140 I = 1, N
               D = D + RS(I)
  140       CONTINUE
            D = D/DBLE(N)
C
C           STEP6  -  UPDATE THETA
C
            THETA = THETA + D
C
C           STEP7  -  STOP IF DESIRED PRECISION IS REACHED
C
            TOL1 = TOL*MAX(1.D0,SIGMB)
            IF (ABS(D).LT.TOL1) THEN
               DSIG = SIGMB - SIGMA
               IF (ABS(DSIG).LT.TOL1) GO TO 160
            END IF
C
C           ITERATIONS CONTINUE IF NIT. LT. MAXIT
C
            SIGMA = SIGMB
            IF (NIT.LT.MAXIT) THEN
               NIT = NIT + 1
               GO TO 60
            ELSE
               IERROR = 5
               NREC = 1
               WRITE (REC,FMT=99990) MAXIT
            END IF
C
C           ITERATIONS ARE TERMINATED
C
  160       CONTINUE
C
C           TEST THAT THE ITERATIONS WERE NOT TERMINATED
C           DUE TO ZERO RESIDUALS
C
            ALLZER = .TRUE.
            DO 180 I = 1, N
               IF (RS(I).NE.ZERO) ALLZER = .FALSE.
  180       CONTINUE
            IF (ALLZER) THEN
               IERROR = 6
               WRITE (REC,FMT=99991)
            END IF
            DO 200 I = 1, N
               RS(I) = X(I) - THETA
  200       CONTINUE
            SIGMA = SIGMB
            GO TO 240
  220       IERROR = 4
            WRITE (REC,FMT=99993) S2
         END IF
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
99995 FORMAT (1X,'** On entry, IPSI.lt.0 or .gt.4: IPSI = ',I16)
99994 FORMAT (1X,'** On entry C.le.0.0 and IPSI=1: C = ',1P,D13.5)
99993 FORMAT (1X,'** Current estimate of SIGMA is zero or negative: SI',
     *       'GMA=',1P,D13.5)
99992 FORMAT (1X,'** All elements of X are equal.')
99991 FORMAT (1X,'** All winsorized residuals are zero.')
99990 FORMAT (1X,'** Number of iterations required exceeds MAXIT: MAXI',
     *       'T =',I12)
99989 FORMAT (1X,'** On entry, IPSI=2 and H1, H2 or H3 has an invalid ',
     *       'value.',/4X,'H1 = ',1P,D13.5,', H2 = ',D13.5,' H3 = ',
     *       D13.5)
99988 FORMAT (1X,'** On entry, IPSI.gt.0 and DCHI.le.0: DCHI = ',1P,
     *       D13.5)
      END
