      SUBROUTINE G02GAF(LINK,MEAN,OFFSET,WEIGHT,N,X,LDX,M,ISX,IP,Y,WT,S,
     *                  A,RSS,IDF,B,IRANK,SE,COV,V,LDV,TOL,MAXIT,IPRINT,
     *                  EPS,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     GLM FOR NORMAL DISTN
C
C
C     V(I,1) = ETA, LINEAR PREDICTOR
C     V(I,2) = FV, FITTED VALUES
C     V(I,3) = VAR, ESTIMATED VARIANCE OF Y(I)
C     V(I,4) = WWT, FINAL VALUE OF WORKING WEIGHT
C     V(I,5) = R, STANDARDIZED RESIDUALS (PRIOR Y)
C     V(I,6) = H, LEVERAGE
C     V(I,7) = OFFSET
C     V(I,8+) = R (IF QR USED) P' (IF SVD USED)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02GAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, EPS, RSS, S, TOL
      INTEGER           IDF, IFAIL, IP, IPRINT, IRANK, LDV, LDX, M,
     *                  MAXIT, N
      CHARACTER         LINK, MEAN, OFFSET, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV(IP*(IP+1)/2), SE(IP), V(LDV,IP+7),
     *                  WK((IP*IP+3*IP+22)/2), WT(*), X(LDX,M), Y(N)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, EPSA, TOLA
      INTEGER           I, IERROR, IFAULT, ITER, J, K, MAXITA, NO, NREC
      LOGICAL           SVD
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G02GAV, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G02GAV, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FBF, G02GAW, G02GBS, G02GBT, G02GBZ, G02GCU,
     *                  G02GCX, G02GCY, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      SVD = .FALSE.
      NREC = 1
      IERROR = 1
      IF (N.LT.2) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (M.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) M
      ELSE IF (IP.LT.1) THEN
         WRITE (P01REC(1),FMT=99978) IP
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99997) LDX, N
      ELSE IF (LDV.LT.N) THEN
         WRITE (P01REC(1),FMT=99988) LDV, N
      ELSE IF (LINK.NE.'I' .AND. LINK.NE.'i' .AND. LINK.NE.'L' .AND.
     *         LINK.NE.'l' .AND. LINK.NE.'R' .AND. LINK.NE.'r' .AND.
     *         LINK.NE.'S' .AND. LINK.NE.'s' .AND. LINK.NE.'E' .AND.
     *         LINK.NE.'e') THEN
         WRITE (P01REC(1),FMT=99989) LINK
      ELSE IF (S.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99993) S
      ELSE IF ((LINK.EQ.'E' .OR. LINK.EQ.'e') .AND. A.EQ.0.0D0) THEN
         WRITE (P01REC(1),FMT=99990)
      ELSE IF (MEAN.NE.'M' .AND. MEAN.NE.'Z' .AND. MEAN.NE.'m' .AND.
     *         MEAN.NE.'z') THEN
         WRITE (P01REC(1),FMT=99996) MEAN
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'U' .AND. WEIGHT.NE.
     *         'w' .AND. WEIGHT.NE.'u') THEN
         WRITE (P01REC(1),FMT=99995) WEIGHT
      ELSE IF (OFFSET.NE.'N' .AND. OFFSET.NE.'Y' .AND. OFFSET.NE.
     *         'n' .AND. OFFSET.NE.'y') THEN
         WRITE (P01REC(1),FMT=99987) OFFSET
      ELSE IF (MAXIT.LT.0) THEN
         WRITE (P01REC(1),FMT=99986) MAXIT
      ELSE IF (TOL.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99985) TOL
      ELSE IF (EPS.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99984) EPS
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
C
C        SET DEFAULT PARAMETERS
C
         ACC = X02AJF()
         IF (MAXIT.EQ.0) THEN
            MAXITA = 10
         ELSE
            MAXITA = MAXIT
         END IF
         IF (TOL.LT.ACC) THEN
            TOLA = ACC*10.0D0
         ELSE
            TOLA = TOL
         END IF
         IF (EPS.LT.ACC) THEN
            EPSA = ACC
         ELSE
            EPSA = EPS
         END IF
C
C        CHECK WEIGHTS
C
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            NO = 0
            DO 20 I = 1, N
               IF (WT(I).LT.0.0D0) GO TO 40
               IF (WT(I).GT.0.0D0) NO = NO + 1
   20       CONTINUE
            GO TO 60
   40       IERROR = 2
            WRITE (P01REC(1),FMT=99994) I
            GO TO 160
         ELSE
            NO = N
         END IF
   60    CONTINUE
         K = 0
         DO 80 J = 1, M
            IF (ISX(J).LT.0) THEN
               IERROR = 3
               WRITE (P01REC(1),FMT=99977) J
               GO TO 160
            ELSE IF (ISX(J).GT.0) THEN
               K = K + 1
            END IF
   80    CONTINUE
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') K = K + 1
         IF (IP.NE.K) THEN
            IERROR = 3
            NREC = 2
            WRITE (P01REC,FMT=99991) IP
            GO TO 160
         ELSE IF (IP.GT.NO) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99992)
            GO TO 160
         END IF
         IF (OFFSET.EQ.'N' .OR. OFFSET.EQ.'n') CALL F06FBF(N,0.0D0,
     *       V(1,7),1)
         CALL F06FBF(N,0.0D0,V(1,6),1)
         CALL F06FBF(N,1.0D0,V(1,3),1)
         CALL G02GCU(N,LINK,A,Y,V(1,2),V(1,1),WT,NO,TOLA,0.0D0)
         CALL G02GBZ(G02GCY,G02GCX,G02GAV,G02GAW,LINK,MEAN,WEIGHT,N,X,
     *               LDX,M,ISX,Y,V(1,6),WT,NO,RSS,IRANK,B,IP,V(1,2),
     *               V(1,1),V(1,3),V(1,4),V(1,7),V(1,8),LDV,A,SVD,TOLA,
     *               MAXITA,ITER,IPRINT,EPSA,WK,IFAULT)
         IF (IFAULT.EQ.1) THEN
            IERROR = 4
            WRITE (P01REC(1),FMT=99982)
            GO TO 160
         ELSE IF (IFAULT.EQ.2) THEN
            IERROR = 5
            WRITE (P01REC(1),FMT=99981)
            GO TO 160
         ELSE IF (IFAULT.EQ.3) THEN
            IERROR = 6
            WRITE (P01REC(1),FMT=99980) MAXITA
         ELSE IF (IFAULT.EQ.4) THEN
            IERROR = 7
            WRITE (P01REC(1),FMT=99983)
         END IF
         IDF = NO - IRANK
         IF (IDF.EQ.0) THEN
            IERROR = 8
            WRITE (P01REC(1),FMT=99979)
         ELSE
            CALL G02GBT(MEAN,N,M,X,LDX,ISX,IP,V(1,8),LDV,SVD,IRANK,
     *                  V(1,4),V(1,6),WK)
            IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
               DO 100 I = 1, N
                  IF (V(I,2).GT.0.0D0) THEN
                     V(I,5) = Y(I) - V(I,2)
                  ELSE
                     V(I,5) = 0.0D0
                  END IF
  100          CONTINUE
            ELSE
               DO 120 I = 1, N
                  IF (WT(I).GT.0.0D0) THEN
                     IF (V(I,2).GT.0.0D0) THEN
                        V(I,5) = SQRT(WT(I))*(Y(I)-V(I,2))
                     ELSE
                        V(I,5) = 0.0D0
                     END IF
                  ELSE
                     V(I,5) = 0.0D0
                  END IF
  120          CONTINUE
            END IF
            CALL G02GBS(IP,IRANK,SVD,V(1,8),LDV,COV,WK)
            IF (S.EQ.0.0D0) S = RSS/IDF
            CALL DSCAL((IP*IP+IP)/2,S,COV,1)
            DO 140 J = 1, IP
               IF (COV((J*J+J)/2).GT.0.0D0) THEN
                  SE(J) = SQRT(COV((J*J+J)/2))
               ELSE
                  SE(J) = 0.0D0
               END IF
  140       CONTINUE
         END IF
      END IF
  160 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, N.lt.2 : N = ',I16)
99998 FORMAT (' ** On entry, M.lt.1 : M = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N : LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, MEAN is not valid : MEAN = ',A1)
99995 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99994 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
99993 FORMAT (' ** On entry, S.lt.0.0 : S = ',D13.5)
99992 FORMAT (' ** Number of requested x-variables greater than N')
99991 FORMAT (' ** On entry, IP incompatible with number of non-zero v',
     *       'alues of ISX:',/'             IP = ',I16)
99990 FORMAT (' ** On entry, A = 0.0 and LINK = E')
99989 FORMAT (' ** On entry, LINK is not valid : LINK = ',A1)
99988 FORMAT (' ** On entry, LDV.lt.N : LDV = ',I16,' N = ',I16)
99987 FORMAT (' ** On entry, OFFSET is not valid : OFFSET = ',A1)
99986 FORMAT (' ** On entry, MAXIT.lt.0 : MAXIT = ',I16)
99985 FORMAT (' ** On entry, TOL.lt.0.0 : TOL = ',D13.5)
99984 FORMAT (' ** On entry, EPS.lt.0.0 : EPS = ',D13.5)
99983 FORMAT (' ** Rank of model has changed during iteration')
99982 FORMAT (' ** Fitted value is at a boundary')
99981 FORMAT (' ** SVD solution failed to converge')
99980 FORMAT (' ** IWLS failed to converge in MAXIT iterations : MAXIT',
     *       ' = ',I16)
99979 FORMAT (' ** Degrees of freedom for error are 0')
99978 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99977 FORMAT (' ** On entry, ISX(',I16,').lt.0')
      END
