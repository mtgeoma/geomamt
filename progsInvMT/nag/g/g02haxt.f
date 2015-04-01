      SUBROUTINE G02HAX(INDW,IPSI,ISIGMA,N,M,X,IX,Y,BETA,CPSI,H1,H2,H3,
     *                  CUCV,D,WGT,THETA,K,SIGMA,RS,TOL,MAXIT,NITMON,
     *                  NIT,WRK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     PURPOSE
C     -------
C     W-ALGORITHM FOR ROBUST AND BOUNDED INFLUENCE LINEAR REGRESSION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BETA, CPSI, CUCV, D, H1, H2, H3, SIGMA, TOL
      INTEGER           IFAIL, INDW, IPSI, ISIGMA, IX, K, M, MAXIT, N,
     *                  NIT, NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  RS(N), THETA(M), WGT(N), WRK((M+4)*N), X(IX,M),
     *                  Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DS, HRS, P0, RPSIP0, SIGMB, SRS, SSW, STD, TAU,
     *                  TMP, TOL1
      INTEGER           I, IFAIL2, IFLAG, J, JN, MED, N0, N4, NITM, NK,
     *                  NN, NOUT
      LOGICAL           SVD
      CHARACTER*80      REC
C     .. External Functions ..
      DOUBLE PRECISION  G02HAU, G07DBX, G07DBY
      EXTERNAL          G02HAU, G07DBX, G07DBY
C     .. External Subroutines ..
      EXTERNAL          F04JGF, F06FCF, M01CAF, DAXPY, DCOPY, DGEMV,
     *                  DSWAP, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, MOD, DBLE, SQRT
C     .. Executable Statements ..
      N4 = 4*N
      K = M
      TAU = MIN(0.000005D0,TOL)
      IFAIL = 0
      NN = (M+1)*M/2
      SIGMB = SIGMA
      NITM = ABS(NITMON)
      IFLAG = 0
      CALL X04ABF(IFLAG,NOUT)
C
C     PRINT OUT HEADINGS FOR ITERATION MONITORING
C
      IF (NITMON.NE.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAF(NOUT,REC)
         WRITE (REC,FMT=99998)
         CALL X04BAF(NOUT,REC)
      END IF
C
C     STEP 1. SET NIT := 1
C     -------
      P0 = G02HAU(0.0D0,IPSI,CPSI,H1,H2,H3)
      NIT = 1
C
C     STEP 2. COMPUTE RESIDUALS R=Y-X*THETA
C     -------
   20 CALL DCOPY(N,Y,1,RS,1)
      CALL DGEMV('N',N,M,-1.0D0,X,IX,THETA,1,1.0D0,RS,1)
C
C     STEP 3. COMPUTE A NEW VALUE SIGMB FOR SIGMA.
C     -------
      IF (ISIGMA.GT.0) THEN
         TMP = 0.0D0
         IF (INDW.EQ.0) THEN
C
C           HUBER-TYPE
C
            N0 = N
            DO 40 I = 1, N
               SRS = RS(I)/SIGMA
               TMP = G07DBY(SRS,IPSI,D) + TMP
   40       CONTINUE
         ELSE
C
C           SCHWEPPE-TYPE
C
            N0 = 0
            DO 60 I = 1, N
               SSW = WGT(I)*SIGMA
               IF (SSW.GT.0.0D0) THEN
                  N0 = N0 + 1
                  SRS = RS(I)/SSW
                  TMP = G07DBY(SRS,IPSI,D)*WGT(I)*WGT(I) + TMP
               END IF
   60       CONTINUE
         END IF
         NK = N0 - K
         IF (NK.LE.0) THEN
            IFAIL = 12
            RETURN
         END IF
         SIGMB = SQRT(TMP/(BETA*DBLE(NK)))*SIGMA
      END IF
      IF (ISIGMA.LT.0) THEN
         IF (INDW.EQ.0) THEN
C
C           HUBER-TYPE
C
            DO 80 I = 1, N
               WRK(I) = ABS(RS(I))
   80       CONTINUE
            N0 = N
         ELSE
C
C           SCHWEPPE-TYPE
C
            N0 = 0
            DO 120 I = 1, N
               IF (WGT(I).LE.0.0D0) GO TO 100
               N0 = N0 + 1
               WRK(N0) = ABS(RS(I))
  100          CONTINUE
  120       CONTINUE
         END IF
         IFAIL2 = 0
         CALL M01CAF(WRK(1),1,N0,'A',IFAIL2)
C
C        IFAIL2 IS DUMMY NO NEED TO CHECK
C
         MED = (N0+1)/2
         SIGMB = WRK(MED)
         IF (2*MED.EQ.N0) SIGMB = (SIGMB+WRK(MED+1))/2.0D0
         SIGMB = SIGMB/BETA
      END IF
C
C     RETURN IF SIGMB LE 0
C
      IF (SIGMB.LE.0.0D0) THEN
         IFAIL = 1
         SIGMA = SIGMB
         RETURN
      END IF
C
C     STEP 4. COMPUTE WEIGHTS AND APPLY THEM TO X; STORE RESULT IN WRK
C     -------
      RPSIP0 = SQRT(P0)
      IF (INDW.EQ.0) THEN
         DO 140 I = 1, N
            IF (RS(I).EQ.0.D0) THEN
               WRK(I) = RPSIP0
               RS(I) = RPSIP0*Y(I)
            ELSE
               SRS = RS(I)/SIGMB
               HRS = SQRT(G07DBX(SRS,IPSI,CPSI,H1,H2,H3)/SRS)
               WRK(I) = HRS
               RS(I) = HRS*Y(I)
            END IF
  140    CONTINUE
      ELSE
         DO 160 I = 1, N
            IF (WGT(I).LE.0.D0) THEN
               WRK(I) = 0.D0
               RS(I) = 0.D0
            ELSE IF (RS(I).EQ.0.D0) THEN
               WRK(I) = RPSIP0
               RS(I) = RPSIP0*Y(I)
            ELSE
               SRS = RS(I)/(SIGMB*WGT(I))
               HRS = SQRT(G07DBX(SRS,IPSI,CPSI,H1,H2,H3)/SRS)
               WRK(I) = HRS
               RS(I) = HRS*Y(I)
            END IF
  160    CONTINUE
      END IF
      DO 180 J = 1, M
         JN = (J+3)*N + 1
         CALL DCOPY(N,X(1,J),1,WRK(JN),1)
         CALL F06FCF(N,WRK,1,WRK(JN),1)
  180 CONTINUE
C
C     STEP 5. SOLVE FOR DELTA
C     -------
      IFAIL2 = 1
      CALL F04JGF(N,M,WRK(N4+1),N,RS,TAU,SVD,STD,K,WRK(1),N4,IFAIL2)
      IF (IFAIL2.EQ.2) THEN
         IFAIL = 2
         RETURN
      END IF
      IF (SVD) IFAIL = 4
C
C     STEP 6. COMPUTE NEW SOLUTION
C     -------
      CALL DSWAP(M,THETA,1,RS,1)
      CALL DAXPY(M,-1.0D0,THETA,1,RS,1)
C
C     ITERATION MONITORING
C
      IF (NITMON.NE.0) THEN
         IF ((MOD(NIT,NITM).EQ.0) .OR. (NIT.EQ.1)) THEN
            J = 1
            WRITE (REC,FMT=99997) NIT, SIGMB, J, THETA(1), RS(1)
            CALL X04BAF(NOUT,REC)
            DO 200 J = 2, M
               WRITE (REC,FMT=99996) J, THETA(J), RS(J)
               CALL X04BAF(NOUT,REC)
  200       CONTINUE
         END IF
      END IF
C
C     STEP 7. STOP ITERATIONS IF DESIRED PRECISION IS REACHED
C     -------
      TOL1 = TOL*SIGMB
      DS = ABS(SIGMA-SIGMB)
      IF (DS.GT.TOL1) GO TO 240
      DO 220 J = 1, M
         TOL1 = ABS(THETA(J))*TOL
         IF (ABS(RS(J)).GT.TOL1) GO TO 240
  220 CONTINUE
      GO TO 260
  240 SIGMA = SIGMB
      IF (NIT.GE.MAXIT) THEN
         IFAIL = 3
         GO TO 260
      END IF
      NIT = NIT + 1
      GO TO 20
  260 SIGMA = SIGMB
C
C     CALCULATE RESIDUALS
C
      CALL DCOPY(N,Y,1,RS,1)
      CALL DGEMV('N',N,M,-1.0D0,X,IX,THETA,1,1.0D0,RS,1)
      RETURN
C
99999 FORMAT ('                ** ITERATION MONITORING FOR THETA **')
99998 FORMAT (' ITERATION       SIGMA        J       THETA       RS')
99997 FORMAT ('   ',I5,'    ',D13.5,'  ',I3,'  ',D13.5,'  ',D13.5)
99996 FORMAT ('                           ',I3,'  ',D13.5,'  ',D13.5)
      END
