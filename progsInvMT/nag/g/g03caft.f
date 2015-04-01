      SUBROUTINE G03CAF(MATRIX,WEIGHT,N,M,X,LDX,NVAR,ISX,NFAC,WT,E,STAT,
     *                  COM,PSI,RES,FL,LDFL,IOP,IWK,WK,LWK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes maximum likelihood estimates of factor loadings
C     The marginal likelihood is maximized using E04LBF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03CAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDFL, LDX, LWK, M, N, NFAC, NVAR
      CHARACTER         MATRIX, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  COM(NVAR), E(NVAR), FL(LDFL,NFAC), PSI(NVAR),
     *                  RES(NVAR*(NVAR-1)/2), STAT(4), WK(LWK), WT(*),
     *                  X(LDX,M)
      INTEGER           IOP(5), ISX(M), IWK(4*NVAR+2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, CHI, COND, CONST, DETS, DF, EPS, ETA, F,
     *                  SCALE, SRN, STEPMX, T, TEMP, WSUM
      INTEGER           I, IERROR, IFAULT, II, IPRINT, IWTP, J, JJ, K,
     *                  K2SPCE, KSPACE, L, LH, LSPACE, MAXFUN, MINWK,
     *                  NO, NREC, NSPACE
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP1(1), WKSP2(1)
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ, G01ECF, DDOT, DNRM2, X02AJF
      INTEGER           P01ABF
      EXTERNAL          F02WDZ, G01ECF, DDOT, DNRM2, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04LBF, F01QCF, F02WUF, F06FBF, F06FCF, F06FDF,
     *                  G03CAW, G03CAX, G03CAY, G03CAZ, DCOPY, F07FDZ,
     *                  DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         INT, LOG, MAX, DBLE, SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') THEN
         MINWK = (5*NVAR*NVAR+33*NVAR-4)/2
      ELSE
         MINWK = MAX((5*NVAR*NVAR+33*NVAR-4)/2,
     *           N*NVAR+7*NVAR+NVAR*(NVAR-1)/2)
      END IF
C
C     check input values
C
      IF (NVAR.LT.2) THEN
         WRITE (P01REC(1),FMT=99999) NVAR
      ELSE IF (N.LE.NVAR) THEN
         WRITE (P01REC(1),FMT=99998) N, NVAR
      ELSE IF (M.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99991) M, NVAR
      ELSE IF ((MATRIX.EQ.'D' .OR. MATRIX.EQ.'S' .OR. MATRIX.EQ.'d' .OR.
     *         MATRIX.EQ.'s') .AND. LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99997) LDX, N
      ELSE IF ((MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') .AND. LDX.LT.M) THEN
         WRITE (P01REC(1),FMT=99992) LDX, M
      ELSE IF (LDFL.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99995) LDFL, NVAR
      ELSE IF (NFAC.LT.1) THEN
         WRITE (P01REC(1),FMT=99988) NFAC
      ELSE IF (NFAC.GT.NVAR) THEN
         WRITE (P01REC(1),FMT=99987) NFAC, NVAR
      ELSE IF (LWK.LT.MINWK) THEN
         WRITE (P01REC(1),FMT=99996) LWK, MINWK
      ELSE IF (MATRIX.NE.'S' .AND. MATRIX.NE.'C' .AND. MATRIX.NE.
     *         's' .AND. MATRIX.NE.'c' .AND. MATRIX.NE.'D' .AND.
     *         MATRIX.NE.'d') THEN
         WRITE (P01REC(1),FMT=99993) MATRIX
      ELSE IF ((MATRIX.EQ.'D' .OR. MATRIX.EQ.'S' .OR. MATRIX.EQ.'d' .OR.
     *         MATRIX.EQ.'s') .AND. (WEIGHT.NE.'W' .AND. WEIGHT.NE.
     *         'w' .AND. WEIGHT.NE.'U' .AND. WEIGHT.NE.'u')) THEN
         WRITE (P01REC(1),FMT=99994) WEIGHT
      ELSE
         IERROR = 0
      END IF
C
C     select and check options
C
      IF (IOP(1).EQ.0) THEN
         IPRINT = -1
         MAXFUN = 100*NVAR
         EPS = X02AJF()
         ACC = 10.0D0*SQRT(EPS)
      ELSE
         IPRINT = IOP(2)
         IF (IPRINT.EQ.0) IPRINT = -1
         MAXFUN = IOP(3)
         ACC = 10.0D0**(-IOP(4))
         EPS = 10.0D0**(-IOP(5))
         IF (MAXFUN.LT.1) THEN
            NREC = 2
            IERROR = 1
            WRITE (P01REC,FMT=99984) MAXFUN
         ELSE IF (ACC.GE.1.0D0 .OR. ACC.LT.X02AJF()) THEN
            IERROR = 1
            WRITE (P01REC(1),FMT=99983) IOP(4)
         ELSE IF (EPS.GE.1.0D0 .OR. EPS.LT.X02AJF()) THEN
            IERROR = 1
            WRITE (P01REC(1),FMT=99982) IOP(5)
         END IF
      END IF
      IF (IERROR.EQ.0) THEN
         NSPACE = 7*NVAR + NVAR*(NVAR-1)/2
C
C        check NVAR
C
         K = 0
         DO 20 I = 1, M
            IF (ISX(I).GT.0) K = K + 1
   20    CONTINUE
         IF (K.NE.NVAR) THEN
            IERROR = 3
            NREC = 2
            WRITE (P01REC,FMT=99990) K, NVAR
            GO TO 560
         END IF
         KSPACE = NSPACE + K*K
         K2SPCE = KSPACE + K*K
         LSPACE = K2SPCE + 6*K - 2
C
         IF (MATRIX.EQ.'D' .OR. MATRIX.EQ.'d' .OR. MATRIX.EQ.'S' .OR.
     *       MATRIX.EQ.'s') THEN
C
C           IF RAW DATA INPUT
C
C           CHECK WEIGHTS
C
            IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
               IWTP = NSPACE + N*(NVAR-1)
               NO = 0
               WSUM = 0.0D0
               SRN = 0.0D0
               DO 40 I = 1, N
                  IF (WT(I).LT.0.0D0) GO TO 60
                  IF (WT(I).GT.0.0D0) THEN
                     NO = NO + 1
                     T = SQRT(WT(I))
                     WK(IWTP+I) = T
                     WSUM = WSUM + WT(I)
                     SRN = SRN + T
                  ELSE
                     WK(IWTP+I) = 0.0D0
                  END IF
   40          CONTINUE
               IF ((MATRIX.EQ.'S' .OR. MATRIX.EQ.'s')
     *             .AND. WSUM.LE.1.0D0) THEN
                  IERROR = 3
                  WRITE (P01REC(1),FMT=99985)
                  GO TO 560
               END IF
            ELSE
               NO = N
               WSUM = DBLE(N)
            END IF
            GO TO 80
   60       IERROR = 2
            WRITE (P01REC(1),FMT=99977) I
            GO TO 560
   80       CONTINUE
            IF (NVAR.GE.INT(WSUM)) THEN
               IERROR = 3
               WRITE (P01REC(1),FMT=99986)
               GO TO 560
            END IF
            CALL G03CAW(WEIGHT,N,X,LDX,M,ISX,NVAR,WT,WSUM,WK(NSPACE+1),
     *                  N,COM,WK)
            SCALE = 1.0D0/(WSUM-1.0D0)
            DO 100 I = 1, K
               COM(I) = COM(I)*SCALE
  100       CONTINUE
C
C           FIND QR DECOMPOSITION
C
            IFAULT = 0
            CALL F01QCF(N,K,WK(NSPACE+1),N,WK,IFAULT)
C
C           CHECK FOR FULL RANK
C
            COND = F02WDZ(K,WK(NSPACE+1),N,WK)
            IF (COND*EPS.GT.1.0D0) THEN
               IERROR = 4
               WRITE (P01REC(1),FMT=99989)
               GO TO 560
            END IF
            DO 120 I = 1, K
               CALL DCOPY(I,WK(NSPACE+(I-1)*N+1),1,WK(NSPACE+(I-1)*K+1),
     *                    1)
  120       CONTINUE
         ELSE
C
C           IF CORRELATION/COVARIANCE MATRIX INPUT
C
            K = 0
            DO 160 I = 1, M
               IF (ISX(I).GT.0) THEN
                  L = NSPACE + K*NVAR
                  K = K + 1
                  COM(K) = X(I,I)
                  DO 140 J = 1, I - 1
                     IF (ISX(J).GT.0) THEN
                        L = L + 1
                        WK(L) = X(J,I)
                     END IF
  140             CONTINUE
                  L = L + 1
                  WK(L) = X(I,I)
               END IF
  160       CONTINUE
C
C           FIND CHOLESKI FACTOR OF CORRELATION/COVARIANCE MATRIX
C
            CALL F07FDZ('U',NVAR,WK(NSPACE+1),NVAR,IFAULT)
            IF (IFAULT.GT.0) THEN
               IERROR = 4
               WRITE (P01REC(1),FMT=99989)
               GO TO 560
            END IF
            WSUM = DBLE(N)
         END IF
C
C        find initial values
C
         CALL F06FBF(K,0.0D0,PSI,1)
         DO 200 I = 1, K
            CALL F06FBF(I-1,0.0D0,FL,1)
            FL(I,1) = 1.0D0
            CALL DTRSV('U','N','N',I,WK(NSPACE+1),K,FL,1)
            DO 180 J = 1, I
               PSI(J) = PSI(J) + FL(J,1)*FL(J,1)
  180       CONTINUE
  200    CONTINUE
         CONST = 1.0D0 - 0.5D0*DBLE(NFAC)/DBLE(K)
         DO 220 I = 1, K
            PSI(I) = CONST/PSI(I)
  220    CONTINUE
C
C        calculate constant for likelihood
C
         DETS = 0.0D0
         DO 240 I = 1, K
            DETS = DETS + LOG(WK(NSPACE+(I-1)*K+I)*WK(NSPACE+(I-1)*K+I))
  240    CONTINUE
C
C        put information in workspace and set bounds
C
         IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') THEN
            CALL F06FDF(K,EPS,COM,1,WK(LSPACE+K+1),1)
            CALL DCOPY(K,COM,1,WK(K2SPCE+K+1),1)
            CALL DCOPY(K,COM,1,WK(LSPACE+2*K+1),1)
         ELSE
            CALL F06FBF(K,EPS,WK(LSPACE+K+1),1)
            CALL F06FBF(K,1.0D0,WK(K2SPCE+K+1),1)
            CALL F06FBF(K,1.0D0,WK(LSPACE+2*K+1),1)
         END IF
         WK(K2SPCE+2*K+1) = DETS
         IWK(K+3) = NFAC
         IWK(K+4) = NSPACE
C
C        set parameters
C
         LH = MAX(K*(K-1)/2,1)
         ETA = 0.9D0
         STEPMX = DBLE(K)
         IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') STEPMX = STEPMX*DNRM2(K,
     *       COM,1)
C
C        call Newton routine
C
         IFAULT = 1
         CALL E04LBF(K,G03CAX,G03CAY,G03CAZ,IPRINT,MAXFUN,ETA,ACC,
     *               STEPMX,0,WK(LSPACE+K+1),WK(LSPACE+2*K+1),PSI,RES,
     *               LH,WK(LSPACE+3*K+1),IWK,F,WK(LSPACE+1),IWK(K+1),4,
     *               WK,LSPACE,IFAULT)
         IF (IFAULT.EQ.-1) THEN
            IERROR = 4
            WRITE (P01REC(1),FMT=99978)
            GO TO 560
         ELSE IF (IFAULT.EQ.-2) THEN
            IERROR = 5
            WRITE (P01REC(1),FMT=99980)
            GO TO 560
         ELSE IF (IFAULT.EQ.2) THEN
            IERROR = 6
            WRITE (P01REC(1),FMT=99981) MAXFUN
         ELSE IF (IFAULT.GT.2) THEN
            IERROR = 7
            WRITE (P01REC(1),FMT=99979)
         END IF
C
C        CALCULATE FACTOR LOADINGS
C
         DO 260 I = 1, K
            IF (PSI(I).NE.WK(LSPACE-K+I)) GO TO 280
  260    CONTINUE
         CALL DCOPY(K,WK(K2SPCE+1),1,E,1)
         GO TO 340
  280    CONTINUE
         DO 300 I = 1, K
            SCALE = 1.0D0/SQRT(PSI(I))
            CALL F06FDF(I,SCALE,WK(NSPACE+(I-1)*K+1),1,WK(KSPACE+(I-1)
     *                  *K+1),1)
  300    CONTINUE
         IFAULT = -1
         CALL F02WUF(K,WK(KSPACE+1),K,0,WKSP1,1,.FALSE.,WKSP2,1,E,
     *               .TRUE.,WK(K2SPCE+2*K+2),IFAULT)
         IF (IFAULT.GT.0) THEN
            IERROR = 5
            WRITE (P01REC(1),FMT=99980)
            GO TO 560
         END IF
         DO 320 I = 1, K
            E(I) = E(I)*E(I)
  320    CONTINUE
  340    CONTINUE
         DO 360 I = 1, K
            WK(K2SPCE+I) = SQRT(PSI(I))
  360    CONTINUE
         DO 380 I = 1, NFAC
            SCALE = E(I) - 1.0D0
            IF (SCALE.GT.0.0D0) THEN
               SCALE = SQRT(SCALE)
            ELSE
               SCALE = 0.0D0
            END IF
            CALL F06FDF(K,SCALE,WK(KSPACE+I),K,FL(1,I),1)
            CALL F06FCF(K,WK(K2SPCE+1),1,FL(1,I),1)
  380    CONTINUE
C
C        calculate residual correlations
C
         IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') THEN
            L = 0
            JJ = 1
            DO 420 J = 2, M
               IF (ISX(J).GT.0) THEN
                  JJ = JJ + 1
                  II = 0
                  DO 400 I = 1, M
                     IF (ISX(I).GT.0) THEN
                        II = II + 1
                        L = L + 1
                        RES(L) = X(I,J) - DDOT(NFAC,FL(II,1),LDFL,
     *                           FL(JJ,1),LDFL)
                        IF (II.GE.JJ-1) GO TO 420
                     END IF
  400             CONTINUE
               END IF
  420       CONTINUE
         ELSE
            L = 0
            DO 460 J = 2, K
               DO 440 I = 1, J - 1
                  L = L + 1
                  RES(L) = DDOT(I,WK(NSPACE+(J-1)*K+1),1,WK(NSPACE+(I-1)
     *                     *K+1),1) - DDOT(NFAC,FL(I,1),LDFL,FL(J,1),
     *                     LDFL)
  440          CONTINUE
  460       CONTINUE
         END IF
C
C        calculate communalities
C
         IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') THEN
            DO 480 I = 1, K
               COM(I) = 1.0D0 - PSI(I)/COM(I)
  480       CONTINUE
         ELSE IF (MATRIX.EQ.'D' .OR. MATRIX.EQ.'d') THEN
            DO 500 I = 1, K
               COM(I) = 1.0D0 - PSI(I)
  500       CONTINUE
C
C           also scale if covariance matrix is required
C
         ELSE
            DO 520 I = 1, K
               WK(K2SPCE+I) = SQRT(COM(I))
               TEMP = 1.0D0 - PSI(I)
               PSI(I) = PSI(I)*COM(I)
               COM(I) = TEMP
  520       CONTINUE
            DO 540 I = 1, NFAC
               CALL F06FCF(K,WK(K2SPCE+1),1,FL(1,I),1)
  540       CONTINUE
         END IF
C
C        calculate test statistics
C
         STAT(1) = F
         DF = 0.5D0*DBLE((NVAR-NFAC)**2-(NVAR+NFAC))
         IF (DF.GT.0.0D0) THEN
            CHI = (WSUM-1.0D0-DBLE(2*NVAR+5)/6.0D0-DBLE(2*NFAC)/3.0D0)*F
            IF (CHI.GE.0.0D0) THEN
               STAT(2) = CHI
               STAT(3) = DF
               IFAULT = 1
               STAT(4) = G01ECF('U',CHI,DF,IFAULT)
            ELSE
               STAT(2) = 0.0D0
               STAT(3) = 0.0D0
               STAT(4) = 0.0D0
            END IF
         ELSE
            STAT(2) = 0.0D0
            STAT(3) = 0.0D0
            STAT(4) = 0.0D0
         END IF
      END IF
  560 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, NVAR.le.1: NVAR = ',I16)
99998 FORMAT (' ** On entry, N.le.NVAR: N = ',I16,' NVAR = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N: LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, LWK is too small: LWK = ',I16,' MIN = ',
     *       I16)
99995 FORMAT (' ** On entry, LDFL.lt.NVAR: LDFL = ',I16,' NVAR = ',I16)
99994 FORMAT (' ** On entry, WEIGHT is not a valid character: WEIGHT = '
     *       ,A1)
99993 FORMAT (' ** On entry, MATRIX is not a valid character: MATRIX = '
     *       ,A1)
99992 FORMAT (' ** On entry, LDX.lt.M: LDX = ',I16,' M = ',I16)
99991 FORMAT (' ** On entry, M.lt.NVAR: M = ',I16,' NVAR = ',I16)
99990 FORMAT (' ** On entry, ',I16,' values of ISX.gt.0, not NVAR:',
     *       /'              NVAR = ',I16)
99989 FORMAT (' ** Matrix not of full rank')
99988 FORMAT (' ** On entry, NFAC.lt.1: NFAC = ',I16)
99987 FORMAT (' ** On entry, NFAC.gt.NVAR: NFAC = ',I16,' NVAR = ',I16)
99986 FORMAT (' ** The number of variables .ge. number of included obs',
     *       'ervations')
99985 FORMAT (' ** The effective number of observations .le. 1')
99984 FORMAT (' ** Maximum number of function evaluations (IOP(3)).le.',
     *       '0:',/'            IOP(3) = ',I16)
99983 FORMAT (' ** Invalid accuracy (IOP(4)) requested: IOP(4) = ',I16)
99982 FORMAT (' ** Invalid eps (IOP(5)) requested: IOP(5) = ',I16)
99981 FORMAT (' ** Estimation has failed to converge in ',I16,' iterat',
     *       'ions')
99980 FORMAT (' ** SVD has failed to converge')
99979 FORMAT (' ** Convergence is uncertain')
99978 FORMAT (' ** Two eigenvalues are equal')
99977 FORMAT (' ** On entry, ',I16,' value of WT.lt.0.0')
      END
