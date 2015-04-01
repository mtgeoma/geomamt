      SUBROUTINE G02GBZ(CLINK,CDER,CDEV,CVAR,LINK,MEAN,WEIGHT,N,X,LDX,M,
     *                  ISX,Y,T,WT,NO,DEV,IRANK,B,IP,FV,ETA,VAR,WWT,OFF,
     *                  Q,LDQ,A,SVD,TOL,MAXIT,ITER,IPRINT,EPS,WK,IERROR)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1659 (JUN 1995).
C
C     COMPUTES ITERATIVE WEIGHTED LEAST SQUARES ESTIMATES FOR
C     GENERALIZED LINEAR MODELS.
C
C     CLINK - LINK FUNCTION
C     CDER  - DERIVATIVES OF LINK FUNCTION
C     CDEV  - DEVIANCE FUNCTION
C     CVAR  - VARIANCE FUNCTION FOR DISTRIBUTION
C
C     ERRORS RETURNED:
C     1 - FITTED VALUE AT BOUNDARY
C     2 - SVD CALCULATION FAILED TO CONVERGE
C     3 - WLS FAILED TO CONVERGE
C     4 - CHANGE IN RANK
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, DEV, EPS, TOL
      INTEGER           IERROR, IP, IPRINT, IRANK, ITER, LDQ, LDX, M,
     *                  MAXIT, N, NO
      LOGICAL           SVD
      CHARACTER         LINK, MEAN, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), ETA(N), FV(N), OFF(N), Q(LDQ,IP), T(N),
     *                  VAR(N), WK(6*IP-5), WT(*), WWT(N), X(LDX,M),
     *                  Y(N)
      INTEGER           ISX(M)
C     .. Function Arguments ..
      DOUBLE PRECISION  CDEV
      EXTERNAL          CDEV
C     .. Subroutine Arguments ..
      EXTERNAL          CDER, CLINK, CVAR
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, DEV1, SQWT
      INTEGER           I, IFAULT, IND, INDQY, IRANK1, J, K, NOUT
      LOGICAL           FINAL
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP(1,1)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ, DDOT
      INTEGER           F06KLF
      EXTERNAL          F02WDZ, DDOT, F06KLF
C     .. External Subroutines ..
      EXTERNAL          F01QCF, F01QDF, F02WUF, F06FBF, F06FCF, DCOPY,
     *                  DGEMV, DTRSV, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MOD, SQRT
C     .. Executable Statements ..
      IERROR = 0
      IF (IPRINT.GT.0) CALL X04ABF(0,NOUT)
      IRANK = IP
      SVD = .FALSE.
      FINAL = .FALSE.
      INDQY = 1
      ITER = 0
C
C     CREATE X MATRIX FOR INITAIL ITERATION
C
      IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
         CALL F06FBF(N,1.0D0,Q,1)
         K = 1
      ELSE
         K = 0
      END IF
      DO 20 J = 1, M
         IF (ISX(J).GT.0) THEN
            K = K + 1
            CALL DCOPY(N,X(1,J),1,Q(1,K),1)
         END IF
   20 CONTINUE
   40 ITER = ITER + 1
C
C     CALCULATE WORKING VARIATE AND WEIGHTS
C
      CALL CDER(N,LINK,ETA,T,WWT,A,WT,NO)
      CALL CVAR(N,FV,T,VAR,WT,NO)
      IF ( .NOT. FINAL) THEN
         DO 60 I = 1, N
            FV(I) = ((ETA(I)-OFF(I))*WWT(I)+Y(I)-FV(I))*VAR(I)
   60    CONTINUE
      END IF
      CALL F06FCF(N,VAR,1,WWT,1)
      IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
         IF (FINAL) THEN
            DO 80 I = 1, N
               IF (WT(I).GT.0.0D0) WWT(I) = WWT(I)*SQRT(WT(I))
   80       CONTINUE
         ELSE
            DO 100 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  SQWT = SQRT(WT(I))
                  WWT(I) = WWT(I)*SQWT
                  FV(I) = FV(I)*SQWT
               END IF
  100       CONTINUE
         END IF
      END IF
      DO 120 I = 1, IP
         CALL F06FCF(N,WWT,1,Q(1,I),1)
  120 CONTINUE
C
C     PERFORM QR
C
      IFAULT = -1
      CALL F01QCF(N,IP,Q,LDQ,WK,IFAULT)
      IF ( .NOT. FINAL) THEN
         IFAULT = -1
         CALL F01QDF('T','S',N,IP,Q,LDQ,WK,1,FV,N,WKSP,IFAULT)
      END IF
      COND = F02WDZ(IP,Q,LDQ,WK)
      IF (COND*EPS.GT.1.0D0) THEN
C
C        PERFORM SVD
C
         SVD = .TRUE.
         IFAULT = -1
         CALL F02WUF(IP,Q,LDQ,INDQY,FV,IP,.FALSE.,WKSP,1,WK,.TRUE.,
     *               WK(IP+1),IFAULT)
         IF (IFAULT.LT.0) THEN
            IERROR = 2
         ELSE
            IRANK = F06KLF(IP,WK(1),1,EPS)
            IF (ITER.EQ.1) THEN
               IRANK1 = IRANK
            ELSE IF (IRANK.NE.IRANK1) THEN
               GO TO 300
C
            END IF
            DO 140 I = 1, IRANK
               WK(I) = 1.0D0/WK(I)
  140       CONTINUE
            DO 160 I = 1, IP
               CALL F06FCF(IRANK,WK,1,Q(1,I),1)
  160       CONTINUE
         END IF
      END IF
      IF (FINAL) THEN
         RETURN
      ELSE
         IF (SVD) THEN
            CALL DGEMV('T',IRANK,IP,1.0D0,Q,LDQ,FV,1,0.0D0,B,1)
         ELSE
            CALL DTRSV('U','N','N',IP,Q,LDQ,FV,1)
            CALL DCOPY(IP,FV,1,B,1)
         END IF
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            K = 1
            CALL F06FBF(N,1.0D0,Q,1)
         ELSE
            K = 0
         END IF
         DO 180 J = 1, M
            IF (ISX(J).GT.0) THEN
               K = K + 1
               CALL DCOPY(N,X(1,J),1,Q(1,K),1)
            END IF
  180    CONTINUE
         IF (N.EQ.NO) THEN
            DO 200 I = 1, N
               ETA(I) = DDOT(IP,B,1,Q(I,1),LDQ) + OFF(I)
  200       CONTINUE
         ELSE
            DO 220 I = 1, N
               IF (WT(I).EQ.0.0D0) THEN
                  ETA(I) = 0.0D0
               ELSE
                  ETA(I) = DDOT(IP,B,1,Q(I,1),LDQ) + OFF(I)
               END IF
  220       CONTINUE
         END IF
         CALL CLINK(N,LINK,ETA,FV,T,A,WT,NO,IND)
         IF (IND.EQ.0) THEN
            GO TO 320
C
         ELSE
            DEV = 0.0D0
            IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
               DO 240 I = 1, N
                  IF (WT(I).GT.0.0D0) THEN
                     DEV = DEV + WT(I)*CDEV(FV(I),Y(I),T(I))
                     IF (T(I).LT.0.0D0) GO TO 320
                  END IF
  240          CONTINUE
            ELSE
               DO 260 I = 1, N
                  DEV = DEV + CDEV(FV(I),Y(I),T(I))
                  IF (T(I).LT.0.0D0) GO TO 320
  260          CONTINUE
            END IF
            IF (DEV.LE.0.0D0) THEN
               FINAL = .TRUE.
               INDQY = 0
            ELSE IF (ITER.EQ.1) THEN
               DEV1 = DEV
            ELSE
               IF (ABS(DEV-DEV1).LT.(1.0D0+DEV)*TOL) THEN
                  FINAL = .TRUE.
                  INDQY = 0
               ELSE
                  DEV1 = DEV
               END IF
            END IF
            IF (IPRINT.GT.0) THEN
               IF (MOD(ITER,IPRINT).EQ.0) THEN
                  WRITE (REC,FMT=99999) ITER, DEV
                  CALL X04BAF(NOUT,REC)
                  IF (SVD) THEN
                     WRITE (REC,FMT=99997)
                     CALL X04BAF(NOUT,REC)
                  END IF
                  DO 280 I = 1, IP
                     WRITE (REC,FMT=99998) I, B(I)
                     CALL X04BAF(NOUT,REC)
  280             CONTINUE
               END IF
            END IF
            IF (ITER.LT.MAXIT .OR. FINAL) THEN
               GO TO 40
C
            ELSE
               IERROR = 3
               FINAL = .TRUE.
               INDQY = 0
               GO TO 40
C
            END IF
         END IF
      END IF
  300 IERROR = 4
      RETURN
  320 IERROR = 1
      RETURN
C
99999 FORMAT (' Iteration ',I16,' Deviance = ',D13.5)
99998 FORMAT (' B(',I6,') = ',D13.5)
99997 FORMAT (' WLS equations are singular')
      END
