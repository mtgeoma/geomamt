      SUBROUTINE G02HAF(INDW,IPSI,ISIGMA,INDC,N,M,X,IX,Y,CPSI,H1,H2,H3,
     *                  CUCV,DCHI,THETA,SIGMA,C,IC,RS,WGT,TOL,MAXIT,
     *                  NITMON,WORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     MASTER ROUTINE FOR ROBUST REGRESSION USING ROBETH ROUTINES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02HAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CPSI, CUCV, DCHI, H1, H2, H3, SIGMA, TOL
      INTEGER           IC, IFAIL, INDC, INDW, IPSI, ISIGMA, IX, M,
     *                  MAXIT, N, NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  C(IC,M), RS(N), THETA(M), WGT(N),
     *                  WORK(4*N+M*(N+M)), X(IX,M), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, SQW, XM, XN
      INTEGER           I, I1, IERROR, IFAIL2, J, K, MM, N4, NIT, NIT1,
     *                  NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02HAT, G02HAX, G02HAY, G02HAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
C     SET UP INITIAL CONSTANTS
C
      NREC = 1
      NIT1 = 0
      NIT = 0
      XM = DBLE(M)
      IERROR = 1
C
C     CHECK FOR ERRORS IN INPUT PARAMETERS
C
      IF (N.LT.2) THEN
         WRITE (REC(1),FMT=99999) N
      ELSE IF (M.LT.1) THEN
         WRITE (REC(1),FMT=99998) M
      ELSE IF (N.LE.M) THEN
         WRITE (REC(1),FMT=99997) N, M
      ELSE IF (IX.LT.N) THEN
         WRITE (REC(1),FMT=99996) N, IX
      ELSE IF (IC.LT.M) THEN
         WRITE (REC(1),FMT=99995) IC, M
      ELSE
         IERROR = 2
      END IF
      IF (IERROR.EQ.1) GO TO 220
      IF (IPSI.LT.0 .OR. IPSI.GT.4) THEN
         WRITE (REC(1),FMT=99988) IPSI
      ELSE
         IERROR = 3
      END IF
      IF (IERROR.EQ.2) GO TO 220
      IF (SIGMA.LE.0.0D0) THEN
         WRITE (REC(1),FMT=99994) SIGMA
      ELSE IF (IPSI.NE.0 .AND. ISIGMA.GT.0 .AND. DCHI.LE.0.0D0) THEN
         WRITE (REC(1),FMT=99993) DCHI
      ELSE IF (IPSI.EQ.1 .AND. CPSI.LE.0.0D0) THEN
         WRITE (REC(1),FMT=99992) CPSI
      ELSE IF (IPSI.EQ.2 .AND. (H1.LT.0.0D0 .OR. H2.LT.H1 .OR. H3.LT.
     *         H2 .OR. (H1.EQ.0.0D0 .AND. H2.EQ.0.0D0 .AND. H3.EQ.0.0D0)
     *         )) THEN
         WRITE (REC(1),FMT=99991)
      ELSE IF (INDW.GT.0 .AND. CUCV.LT.SQRT(XM)) THEN
         WRITE (REC(1),FMT=99990) CUCV
      ELSE IF (INDW.LT.0 .AND. CUCV.LT.XM) THEN
         WRITE (REC(1),FMT=99989) CUCV
      ELSE
         IERROR = 4
      END IF
      IF (IERROR.EQ.3) GO TO 220
      IF (TOL.LE.0.0D0) THEN
         WRITE (REC(1),FMT=99987) TOL
      ELSE IF (MAXIT.LE.0) THEN
         WRITE (REC(1),FMT=99986) MAXIT
      ELSE
         IERROR = 12
      END IF
      IF (IERROR.EQ.4) GO TO 220
C
C      SET UP VALUES OF CONSTANTS
C
      XN = DBLE(N)
      MM = (M+1)*M/2
      N4 = N*4
C
C     CALCULATE WEIGHTS
C
      IF (INDW.NE.0) THEN
         IFAIL2 = 0
         CALL G02HAY(INDW,N,M,X,IX,MM,CUCV,WGT,MAXIT,NITMON,TOL,NIT1,
     *               WORK(1),WORK(MM+1),WORK(MM+N+1),IFAIL2)
         IF (IFAIL2.EQ.1) THEN
            IERROR = 5
            NREC = 2
            WRITE (REC,FMT=99985) MAXIT
            GO TO 220
         END IF
      ELSE
         DO 20 I = 1, N
            WGT(I) = 1.0D0
   20    CONTINUE
      END IF
C
C
C     CALCULATE VALUE OF BETA
C
      IFAIL2 = 0
      CALL G02HAT(ISIGMA,INDW,IPSI,N,DCHI,WGT,BETA,MAXIT,TOL,WORK,
     *            IFAIL2)
      IF (IFAIL2.EQ.1) THEN
         IERROR = 6
         NREC = 2
         WRITE (REC,FMT=99984) MAXIT
      END IF
C
C     TRANSFORM FOR MALLOWS TYPE ESTIMATOR
C
      IF (INDW.LT.0) THEN
         DO 60 I = 1, N
            SQW = SQRT(WGT(I))
            WGT(I) = SQW
            Y(I) = Y(I)*SQW
            DO 40 J = 1, M
               X(I,J) = X(I,J)*SQW
   40       CONTINUE
   60    CONTINUE
      END IF
C
C     ESTIMATE PARAMETERS USING ITERATIVE-WEIGHTED ALGORITHM
C
      IFAIL2 = 0
      CALL G02HAX(INDW,IPSI,ISIGMA,N,M,X,IX,Y,BETA,CPSI,H1,H2,H3,CUCV,
     *            DCHI,WGT,THETA,K,SIGMA,RS,TOL,MAXIT,NITMON,NIT,WORK,
     *            IFAIL2)
C
C     CALCULATE VARIANCE OF ESTIMATES
C
C     BACK TRANSFORM FOR MALLOWS TYPE ESTIMATORS
C
      IF (INDW.LT.0) THEN
         DO 100 J = 1, M
            DO 80 I = 1, N
               IF (WGT(I).GT.0.0D0) THEN
                  X(I,J) = X(I,J)/WGT(I)
               END IF
   80       CONTINUE
  100    CONTINUE
         DO 120 I = 1, N
            IF (WGT(I).GT.0.0D0) THEN
               RS(I) = RS(I)/WGT(I)
               Y(I) = Y(I)/WGT(I)
               WGT(I) = WGT(I)*WGT(I)
            END IF
  120    CONTINUE
      END IF
      IF (IFAIL2.EQ.12) THEN
         WRITE (REC,FMT=99976) N, K
      ELSE IF (IFAIL2.EQ.1) THEN
         WRITE (REC(1),FMT=99983)
      ELSE
         IERROR = 7
      END IF
      IF (IERROR.EQ.12) GO TO 200
      IF (IFAIL2.EQ.2) THEN
         WRITE (REC(1),FMT=99982)
      ELSE IF (IFAIL2.EQ.3) THEN
         NREC = 2
         WRITE (REC,FMT=99981) MAXIT
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.7) GO TO 200
      IF (IFAIL2.EQ.4) THEN
         WRITE (REC(1),FMT=99980) K
         IERROR = 8
         GO TO 200
      END IF
      IFAIL2 = 0
      CALL G02HAZ(INDW,IPSI,INDC,SIGMA,N,M,X,IX,RS,WGT,CPSI,H1,H2,H3,C,
     *            IC,WORK(1),WORK(N+1),WORK(2*N+1),IFAIL2)
      IF (IFAIL2.EQ.1) THEN
         WRITE (REC(1),FMT=99979)
         IERROR = 9
         GO TO 200
      END IF
      IF (IFAIL2.EQ.2) THEN
         WRITE (REC(1),FMT=99978)
         IERROR = 10
         GO TO 200
      END IF
C
C     CALCULATE CORRELATIONS AND SE
C
      IERROR = 0
      DO 140 I = 1, M
         WORK(I) = 1
         IF (C(I,I).LE.0.0D0) THEN
            IERROR = 11
            WRITE (REC(1),FMT=99977)
            WORK(I) = -1.0D0
         ELSE
            C(I,I) = SQRT(C(I,I))
         END IF
  140 CONTINUE
      DO 180 I = 2, M
         I1 = I - 1
         DO 160 J = 1, I1
            IF (WORK(I).GT.0.0D0 .AND. WORK(J).GT.0.0D0) THEN
               C(J,I) = C(J,I)/(C(I,I)*C(J,J))
            ELSE
               C(J,I) = 0.0D0
            END IF
  160    CONTINUE
  180 CONTINUE
C
C     ADD INFORMATION TO WORK ARRAY
C
  200 WORK(1) = BETA
      WORK(2) = DBLE(NIT1)
      WORK(3) = DBLE(NIT)
      WORK(4) = DBLE(K)
  220 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** On entry, N lt 2: N=',I16)
99998 FORMAT (' ** On entry, M lt 1: M=',I16)
99997 FORMAT (' ** On entry, N le M: N=',I16,' M=',I16)
99996 FORMAT (' ** On entry, IX lt N: N=',I16,' IX=',I16)
99995 FORMAT (' ** On entry, IC lt M: IC=',I16,'M=',I16)
99994 FORMAT (' ** On entry, SIGMA le 0: SIGMA=',D13.5)
99993 FORMAT (' ** On entry, DCHI le 0: DCHI=',D13.5)
99992 FORMAT (' ** On entry, CPSI value not permitted: CPSI=',D13.5)
99991 FORMAT (' ** On entry, H1,H2,H3 incorrectly set')
99990 FORMAT (' ** On entry, CUCV lt SQRT(M): CUCV=',D13.5)
99989 FORMAT (' ** On entry, CUCV lt M: CUCV=',D13.5)
99988 FORMAT (' ** On entry, IPSI does not have a valid value: IPSI=',
     *       I16)
99987 FORMAT (' ** On entry, TOL le 0: TOL=',D13.5)
99986 FORMAT (' ** On entry, MAXIT le 0: MAXIT=',I16)
99985 FORMAT (' ** Iterations to calculate weights failed to converge',
     *       /'    in MAXIT iterations: MAXIT = ',I16)
99984 FORMAT (' ** Iterations to calculate estimates of BETA failed to',
     *       ' converge  ',/'    in MAXIT iterations: MAXIT = ',I16)
99983 FORMAT (' ** Estimated value of SIGMA is zero')
99982 FORMAT (' ** Iterations to solve weighted least squares estimate',
     *       's   failed to converge')
99981 FORMAT (' ** Iterations to calculate estimates of THETA failed t',
     *       'o converge ',/'    in MAXIT iterations: MAXIT = ',I16)
99980 FORMAT (' ** Weighted least squares equations not of full rank, ',
     *       'rank=',I16)
99979 FORMAT (' ** Failure to invert matrix while calculating covarian',
     *       'ce')
99978 FORMAT (' ** Factor for covariance matrix = 0, uncorrected INV(T',
     *       'RAN(X)*X) given ')
99977 FORMAT (' ** Variance of an element of THETA le 0, correlations ',
     *       'set  to 0 ')
99976 FORMAT (1X,'** Error degrees of freedom le 0: N = ',I16,' Rank o',
     *       'f X = ',I16)
      END
