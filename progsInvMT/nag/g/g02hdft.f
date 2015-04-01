      SUBROUTINE G02HDF(CHI,PSI,PSIP0,BETA,INDW,ISIGMA,N,M,X,IX,Y,WGT,
     *                  THETA,K,SIGMA,RS,TOL,EPS,MAXIT,NITMON,NIT,WK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     PURPOSE
C     -------
C     IWLS-ALGORITHM FOR ROBUST AND BOUNDED INFLUENCE LINEAR REGRESSION
C     BASED ON ROUTINES IN ROBETH BY A MARAZZI
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02HDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BETA, EPS, PSIP0, SIGMA, TOL
      INTEGER           IFAIL, INDW, ISIGMA, IX, K, M, MAXIT, N, NIT,
     *                  NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  RS(N), THETA(M), WGT(N), WK((M+4)*N), X(IX,M),
     *                  Y(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  CHI, PSI
      EXTERNAL          CHI, PSI
C     .. Local Scalars ..
      DOUBLE PRECISION  HRS, RPSIP0, SIGMB, SQW, SRS, SSW, STD, TMP,
     *                  TOL1, U
      INTEGER           I, IERROR, IFAIL2, IFLAG, J, JN, MED, N0, N4,
     *                  NK, NOUT, NREC
      LOGICAL           ITRM, SVD
      CHARACTER*80      REC
C     .. Local Arrays ..
      CHARACTER*80      EREC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04JGF, F06FCF, M01CAF, DAXPY, DCOPY, DGEMV,
     *                  DSWAP, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MOD, DBLE, SQRT
C     .. Executable Statements ..
C
C     CHECK FOR ERRORS IN INPUT PARAMETERS
C
      NIT = 0
      NREC = 1
      ITRM = .FALSE.
      IERROR = 1
      IF (N.LT.2) THEN
         WRITE (EREC(1),FMT=99995) N
      ELSE IF (M.LT.1) THEN
         WRITE (EREC(1),FMT=99994) M
      ELSE IF (N.LE.M) THEN
         WRITE (EREC(1),FMT=99993) N, M
      ELSE IF (IX.LT.N) THEN
         WRITE (EREC(1),FMT=99992) N, IX
      ELSE
         IERROR = 2
      END IF
      IF (IERROR.NE.1) THEN
         IF (SIGMA.LE.0.0D0) THEN
            WRITE (EREC(1),FMT=99991) SIGMA
         ELSE
            IERROR = 3
         END IF
         IF (ISIGMA.NE.0 .AND. IERROR.EQ.3) THEN
            IF (BETA.LE.0.0D0) THEN
               WRITE (EREC(1),FMT=99990) BETA
               IERROR = 2
            END IF
         END IF
         IF (IERROR.NE.2) THEN
            IF (TOL.LE.0.0D0) THEN
               WRITE (EREC(1),FMT=99989) TOL
            ELSE IF (MAXIT.LE.0) THEN
               WRITE (EREC(1),FMT=99988) MAXIT
            ELSE
               IERROR = 0
            END IF
            IF (IERROR.NE.3) THEN
C
               K = M
               N4 = 4*N
               SIGMB = SIGMA
C
C              PRINT OUT HEADINGS FOR ITERATION MONITORING
C
               IF (NITMON.GT.0) THEN
                  IFLAG = 0
                  CALL X04ABF(IFLAG,NOUT)
                  WRITE (REC,FMT=99999)
                  CALL X04BAF(NOUT,REC)
                  WRITE (REC,FMT=99998)
                  CALL X04BAF(NOUT,REC)
               END IF
C
C              TRANSFORM FOR MALLOWS TYPE REGRESSION
C
               IF (INDW.LT.0) THEN
                  ITRM = .TRUE.
                  DO 40 I = 1, N
                     IF (WGT(I).GT.0.0D0) THEN
                        SQW = SQRT(WGT(I))
                        WGT(I) = SQW
                        Y(I) = Y(I)*SQW
                        DO 20 J = 1, M
                           X(I,J) = X(I,J)*SQW
   20                   CONTINUE
                     END IF
   40             CONTINUE
               END IF
C
C
C              STEP 1. SET NIT := 1
C              -------
               NIT = 1
   60          CONTINUE
C
C              STEP 2. COMPUTE RESIDUALS R=Y-X*THETA
C              -------
C
               CALL DCOPY(N,Y,1,RS,1)
               CALL DGEMV('N',N,M,-1.0D0,X,IX,THETA,1,1.0D0,RS,1)
C
C
C              STEP 3. COMPUTE A NEW VALUE SIGMB FOR SIGMA.
C              -------
               IF (ISIGMA.GT.0) THEN
                  TMP = 0.0D0
                  IF (INDW.EQ.0) THEN
C
C                    HUBER-TYPE
C
                     N0 = N
                     DO 80 I = 1, N
                        SRS = RS(I)/SIGMA
                        U = CHI(SRS)
                        IF (U.LT.0.0D0) THEN
                           GO TO 420
C
                        ELSE
C
                           TMP = U + TMP
                        END IF
   80                CONTINUE
                  ELSE
C
C                    SCHWEPPE-TYPE
C
                     N0 = 0
                     DO 100 I = 1, N
                        SSW = WGT(I)*SIGMA
                        IF (SSW.GT.0.0D0) THEN
                           N0 = N0 + 1
                           SRS = RS(I)/SSW
                           U = CHI(SRS)
                           IF (U.LT.0.0D0) THEN
                              GO TO 280
C
                           ELSE
C
                              TMP = WGT(I)*U*WGT(I) + TMP
                           END IF
                        END IF
  100                CONTINUE
                  END IF
                  NK = N0 - K
                  IF (NK.LE.0) THEN
                     GO TO 400
C
                  ELSE
C
                     SIGMB = SQRT(TMP/(DBLE(NK)*BETA))*SIGMA
                  END IF
               END IF
               IF (ISIGMA.LT.0) THEN
                  IF (INDW.EQ.0) THEN
C
C                    USING MEDIAN ABSOLUTE DEVIATION TO ESTIMATE SIGMA
C
C                    HUBER-TYPE
C
                     DO 120 I = 1, N
                        WK(I) = ABS(RS(I))
  120                CONTINUE
                     N0 = N
                  ELSE
C
C                    SCHWEPPE-TYPE
C
                     N0 = 0
                     DO 140 I = 1, N
                        IF (WGT(I).GT.0.0D0) THEN
                           N0 = N0 + 1
                           WK(N0) = ABS(RS(I))
                        END IF
  140                CONTINUE
                  END IF
                  IFAIL2 = 0
                  CALL M01CAF(WK,1,N0,'A',IFAIL2)
C
C                 IFAIL2 IS DUMMY NO NEED TO CHECK
C
                  MED = (N0+1)/2
                  SIGMB = WK(MED)
                  IF (2*MED.EQ.N0) SIGMB = (SIGMB+WK(MED+1))/2.0D0
                  SIGMB = SIGMB/BETA
               END IF
C
C              RETURN IF SIGMB LE 0
C
               IF (SIGMB.LE.0.0D0) THEN
                  GO TO 380
C
               ELSE
C
C
C                 STEP 4. COMPUTE WEIGHTS AND APPLY THEM TO X;
C                 -------
C                         STORE RESULT IN WK
C
                  RPSIP0 = SQRT(PSIP0)
                  IF (INDW.EQ.0) THEN
                     DO 160 I = 1, N
                        IF (RS(I).EQ.0.0D0) THEN
                           HRS = RPSIP0
                        ELSE
                           SRS = RS(I)/SIGMB
                           HRS = SQRT(PSI(SRS)/SRS)
                        END IF
                        WK(I) = HRS
                        RS(I) = Y(I)*HRS
  160                CONTINUE
                  ELSE
                     DO 180 I = 1, N
                        IF (WGT(I).LE.0.0D0) THEN
                           HRS = 0.0D0
                        ELSE IF (RS(I).EQ.0.0D0) THEN
                           HRS = RPSIP0
                        ELSE
                           SRS = RS(I)/(WGT(I)*SIGMB)
                           HRS = SQRT(PSI(SRS)/SRS)
                        END IF
                        WK(I) = HRS
                        RS(I) = Y(I)*HRS
  180                CONTINUE
                  END IF
                  DO 200 J = 1, M
                     JN = (J+3)*N + 1
                     CALL DCOPY(N,X(1,J),1,WK(JN),1)
                     CALL F06FCF(N,WK,1,WK(JN),1)
  200             CONTINUE
C
C                 STEP 5. SOLVE FOR DELTA
C                 -------
C
                  IFAIL2 = 1
                  CALL F04JGF(N,M,WK(N4+1),N,RS,EPS,SVD,STD,K,WK,N4,
     *                        IFAIL2)
                  IF (IFAIL2.EQ.2) THEN
                     GO TO 360
C
                  ELSE
C
                     IF (SVD) THEN
                        IERROR = 7
                        WRITE (EREC(1),FMT=99983) K
                     END IF
C
C                    STEP 6. PLACE NEW SOLUTION IN THETA
C                    -------
                     CALL DSWAP(M,THETA,1,RS,1)
                     CALL DAXPY(M,-1.0D0,THETA,1,RS,1)
C
C                    ITERATION MONITORING
C
                     IF (NITMON.GT.0) THEN
                        IF ((MOD(NIT,NITMON).EQ.0) .OR. (NIT.EQ.1)) THEN
                           J = 1
                           WRITE (REC,FMT=99997) NIT, SIGMB, J,
     *                       THETA(1), RS(1)
                           CALL X04BAF(NOUT,REC)
                           DO 220 J = 2, M
                              WRITE (REC,FMT=99996) J, THETA(J), RS(J)
                              CALL X04BAF(NOUT,REC)
  220                      CONTINUE
                        END IF
                     END IF
C
C                    STEP 7. STOP ITERATIONS IF DESIRED PRECISION
C                    -------
C                            IS REACHED
C
                     IF (ISIGMA.NE.0) THEN
                        TOL1 = TOL*SIGMB
                        IF (ABS(SIGMA-SIGMB).GT.TOL1) GO TO 260
                     END IF
                     DO 240 J = 1, M
                        TOL1 = ABS(THETA(J))*TOL
                        IF (ABS(RS(J)).GT.TOL1) GO TO 260
  240                CONTINUE
                     GO TO 300
C
  260                SIGMA = SIGMB
                     IF (NIT.GE.MAXIT) THEN
                        GO TO 320
C
                     ELSE
C
                        NIT = NIT + 1
                        GO TO 60
C
                     END IF
                  END IF
               END IF
  280          IERROR = 4
               WRITE (EREC(1),FMT=99987) SRS, U
               GO TO 440
C
  300          SIGMA = SIGMB
               GO TO 340
C
  320          IERROR = 8
               NREC = 2
               WRITE (EREC,FMT=99984) MAXIT
C
C              CALCULATE RESIDUALS
C
  340          CALL DCOPY(N,Y,1,RS,1)
               CALL DGEMV('N',N,M,-1.0D0,X,IX,THETA,1,1.0D0,RS,1)
               GO TO 440
C
  360          IERROR = 6
               WRITE (EREC(1),FMT=99985)
               GO TO 440
C
  380          IERROR = 5
               WRITE (EREC(1),FMT=99986)
               GO TO 440
C
  400          IERROR = 9
               WRITE (EREC(1),FMT=99982) K, N0
               GO TO 440
C
  420          IERROR = 4
               WRITE (EREC(1),FMT=99987) SRS, U
            END IF
         END IF
      END IF
C
C     BACK TRANSFORM IN MALLOWS CASE
C
  440 IF (ITRM) THEN
         DO 480 I = 1, N
            IF (WGT(I).GT.0.0D0) THEN
               Y(I) = Y(I)/WGT(I)
               RS(I) = RS(I)/WGT(I)
               DO 460 J = 1, M
                  X(I,J) = X(I,J)/WGT(I)
  460          CONTINUE
               WGT(I) = WGT(I)*WGT(I)
            ELSE
               RS(I) = 0.0D0
            END IF
  480    CONTINUE
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,EREC)
      RETURN
C
99999 FORMAT ('                    ** ITERATION MONITORING **')
99998 FORMAT (' ITERATION       SIGMA        J       THETA       DELTA')
99997 FORMAT ('   ',I5,'    ',D13.5,'  ',I3,'  ',D13.5,'  ',D13.5)
99996 FORMAT ('                           ',I3,'  ',D13.5,'  ',D13.5)
99995 FORMAT (' ** On entry, N lt 2: N=',I16)
99994 FORMAT (' ** On entry, M lt 1: M=',I16)
99993 FORMAT (' ** On entry, N le M: N=',I16,' M=',I16)
99992 FORMAT (' ** On entry, IX lt N: N=',I16,' IX=',I16)
99991 FORMAT (' ** On entry, SIGMA le 0: SIGMA=',D13.5)
99990 FORMAT (' ** On entry, BETA le 0: BETA=',D13.5)
99989 FORMAT (' ** On entry, TOL le 0: TOL=',D13.5)
99988 FORMAT (' ** On entry, MAXIT le 0: MAXIT=',I16)
99987 FORMAT (' ** Value given by CHI function lt 0: CHI(',D13.5,') =',
     *       D13.5)
99986 FORMAT (' ** Estimated value of SIGMA is zero')
99985 FORMAT (' ** Iterations to solve weighted least squares equation',
     *       's failed   to converge')
99984 FORMAT (' ** Iterations to calculate estimates of THETA failed t',
     *       'o ',/'    converge in MAXIT iterations: MAXIT = ',I16)
99983 FORMAT (' ** Weighted least squares equations not of full rank: ',
     *       'rank=',I16)
99982 FORMAT (' Value of N-K le 0: N= ',I16,' K= ',I16)
      END
