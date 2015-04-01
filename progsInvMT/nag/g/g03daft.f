      SUBROUTINE G03DAF(WEIGHT,N,M,X,LDX,ISX,NVAR,ING,NG,WT,NIG,GMEAN,
     *                  LDG,DET,GC,STAT,DF,SIG,WK,IWK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes test for equality of within group covariance matrices
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03DAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DF, SIG, STAT
      INTEGER           IFAIL, LDG, LDX, M, N, NG, NVAR
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  DET(NG), GC((NG+1)*NVAR*(NVAR+1)/2),
     *                  GMEAN(LDG,NVAR), WK(N*(NVAR+1)), WT(*), X(LDX,M)
      INTEGER           ING(N), ISX(M), IWK(NG), NIG(NG)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, DETP, RSUM, SUM, TEMP, WSUM
      INTEGER           I, IERROR, IFAULT, J, K, L, NREC, NWT
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01QCF, F06FDF, G03DAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG, DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 1
      IF (NVAR.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) NVAR
      ELSE IF (NG.LT.2) THEN
         WRITE (P01REC(1),FMT=99998) NG
      ELSE IF (N.LT.1) THEN
         WRITE (P01REC(1),FMT=99997) N
      ELSE IF (M.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99996) M, NVAR
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99995) LDX, N
      ELSE IF (LDG.LT.NG) THEN
         WRITE (P01REC(1),FMT=99994) LDG, NG
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u') THEN
         WRITE (P01REC(1),FMT=99993) WEIGHT
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
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
            WRITE (P01REC,FMT=99992) K, NVAR
            GO TO 800
         END IF
         DO 40 I = 1, NG
            NIG(I) = 0
   40    CONTINUE
C
C        check ING (and WT)
C
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            DO 60 I = 1, NG
               DET(I) = 0.0D0
   60       CONTINUE
            WSUM = 0.0D0
            DO 80 I = 1, N
               IF (WT(I).LT.0.0D0) THEN
                  IERROR = 2
                  WRITE (P01REC(1),FMT=99990) I
                  GO TO 800
               ELSE IF (WT(I).NE.0.0D0) THEN
                  K = ING(I)
                  IF (K.LE.0 .OR. K.GT.NG) THEN
                     IERROR = 3
                     WRITE (P01REC(1),FMT=99991) I
                     GO TO 800
                  END IF
                  NIG(K) = NIG(K) + 1
                  DET(K) = DET(K) + WT(I)
                  WSUM = WSUM + WT(I)
               END IF
   80       CONTINUE
         ELSE
            DO 100 I = 1, N
               K = ING(I)
               IF (ING(I).LE.0 .OR. ING(I).GT.NG) THEN
                  IERROR = 3
                  WRITE (P01REC(1),FMT=99991) I
                  GO TO 800
               END IF
               NIG(K) = NIG(K) + 1
  100       CONTINUE
            WSUM = N
            DO 120 I = 1, NG
               DET(I) = DBLE(NIG(I))
  120       CONTINUE
         END IF
         DO 140 I = 1, NG
            IF (DET(I).LT.1.0D0) THEN
               IERROR = 3
               NREC = 2
               WRITE (P01REC,FMT=99989) I
               GO TO 800
            ELSE IF (NIG(I).LT.NVAR) THEN
               IERROR = 3
               NREC = 2
               WRITE (P01REC,FMT=99988) I
               GO TO 800
            END IF
  140    CONTINUE
         NWT = N*NVAR
C
C        move data to workspace
C
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            K = 0
            DO 200 J = 1, M
               IF (ISX(J).GT.0) THEN
                  IWK(1) = 1
                  DO 160 I = 1, NG - 1
                     IWK(I+1) = IWK(I) + NIG(I)
  160             CONTINUE
                  DO 180 I = 1, N
                     IF (WT(I).NE.0.0D0) THEN
                        L = ING(I)
                        WK(N*K+IWK(L)) = X(I,J)
                        IWK(L) = IWK(L) + 1
                     END IF
  180             CONTINUE
                  K = K + 1
               END IF
  200       CONTINUE
            IWK(1) = 1
            DO 220 I = 1, NG - 1
               IWK(I+1) = IWK(I) + NIG(I)
  220       CONTINUE
            DO 240 I = 1, N
               IF (WT(I).NE.0.0D0) THEN
                  L = ING(I)
                  WK(NWT+IWK(L)) = WT(I)
                  IWK(L) = IWK(L) + 1
               END IF
  240       CONTINUE
         ELSE
            K = 0
            DO 300 J = 1, M
               IF (ISX(J).GT.0) THEN
                  IWK(1) = 1
                  DO 260 I = 1, NG - 1
                     IWK(I+1) = IWK(I) + NIG(I)
  260             CONTINUE
                  DO 280 I = 1, N
                     L = ING(I)
                     WK(N*K+IWK(L)) = X(I,J)
                     IWK(L) = IWK(L) + 1
  280             CONTINUE
                  K = K + 1
               END IF
  300       CONTINUE
         END IF
C
C        compute means
C
         IWK(1) = NIG(1)
         DO 320 I = 2, NG
            IWK(I) = IWK(I-1) + NIG(I)
  320    CONTINUE
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            DO 400 J = 1, NVAR
               SUM = 0.0D0
               L = N*(J-1)
               DO 340 I = 1, IWK(1)
                  SUM = SUM + WK(L+I)*WK(NWT+I)
  340          CONTINUE
               GMEAN(1,J) = SUM/DET(1)
               DO 380 K = 2, NG
                  SUM = 0.0D0
                  DO 360 I = IWK(K-1) + 1, IWK(K)
                     SUM = SUM + WK(L+I)*WK(NWT+I)
  360             CONTINUE
                  GMEAN(K,J) = SUM/DET(K)
  380          CONTINUE
  400       CONTINUE
         ELSE
            DO 480 J = 1, NVAR
               SUM = 0.0D0
               L = N*(J-1)
               DO 420 I = 1, IWK(1)
                  SUM = SUM + WK(L+I)
  420          CONTINUE
               GMEAN(1,J) = SUM/DET(1)
               DO 460 K = 2, NG
                  SUM = 0.0D0
                  DO 440 I = IWK(K-1) + 1, IWK(K)
                     SUM = SUM + WK(L+I)
  440             CONTINUE
                  GMEAN(K,J) = SUM/DET(K)
  460          CONTINUE
  480       CONTINUE
         END IF
C
C        mean centre data in workspace
C
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            DO 500 I = 1, IWK(NG)
               IF (WK(NWT+I).GT.0.0D0) WK(NWT+I) = SQRT(WK(NWT+I))
  500       CONTINUE
            DO 580 J = 1, NVAR
               L = N*(J-1)
               DO 520 I = 1, IWK(1)
                  WK(L+I) = WK(NWT+I)*(WK(L+I)-GMEAN(1,J))
  520          CONTINUE
               DO 560 K = 2, NG
                  DO 540 I = IWK(K-1) + 1, IWK(K)
                     WK(L+I) = WK(NWT+I)*(WK(L+I)-GMEAN(K,J))
  540             CONTINUE
  560          CONTINUE
  580       CONTINUE
         ELSE
            DO 660 J = 1, NVAR
               L = N*(J-1)
               DO 600 I = 1, IWK(1)
                  WK(L+I) = WK(L+I) - GMEAN(1,J)
  600          CONTINUE
               DO 640 K = 2, NG
                  DO 620 I = IWK(K-1) + 1, IWK(K)
                     WK(L+I) = WK(L+I) - GMEAN(K,J)
  620             CONTINUE
  640          CONTINUE
  660       CONTINUE
         END IF
C
C        find QR decomposition of the variables for each group
C
         NWT = NWT + 1
         K = 1
         L = NVAR*(NVAR+1)/2 + 1
         DO 700 I = 1, NG
            IFAULT = -1
            CALL F01QCF(NIG(I),NVAR,WK(K),N,WK(NWT),IFAULT)
            TEMP = 1.0D0/SQRT(DET(I)-1.0D0)
            DO 680 J = 1, NVAR
               CALL F06FDF(J,TEMP,WK(N*(J-1)+K),1,GC(L),1)
               L = L + J
  680       CONTINUE
            K = K + NIG(I)
  700    CONTINUE
C
C        find QR decomposition for pooled matrix
C
         CALL G03DAZ(NG,IWK,NVAR,WK,N)
         L = 1
         TEMP = 1.0D0/SQRT(WSUM-DBLE(NG))
         DO 720 J = 1, NVAR
            CALL F06FDF(J,TEMP,WK(N*(J-1)+1),1,GC(L),1)
            L = L + J
  720    CONTINUE
C
C        compute determinants
C
         L = 0
         SUM = 0.0D0
         DO 740 I = 1, NVAR
            L = L + I
            TEMP = GC(L)
            IF (TEMP.EQ.0.0D0) THEN
               IERROR = 4
               WRITE (P01REC,FMT=99986)
               GO TO 800
            END IF
            SUM = SUM + LOG(ABS(TEMP))
  740    CONTINUE
         DETP = 2.0D0*SUM
         WSUM = WSUM - DBLE(NG)
         STAT = WSUM*DETP
         WK(1) = DETP
         RSUM = 0.0D0
         DO 780 J = 1, NG
            SUM = 0.0D0
            DO 760 I = 1, NVAR
               L = L + I
               TEMP = GC(L)
               IF (TEMP.EQ.0.0D0) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99987) J
                  GO TO 800
               END IF
               SUM = SUM + LOG(ABS(TEMP))
  760       CONTINUE
            TEMP = DET(J) - 1.0D0
            RSUM = RSUM + 1.0D0/TEMP
            STAT = STAT - TEMP*2.0D0*SUM
            DET(J) = 2.0D0*SUM
  780    CONTINUE
C
C        compute test statistics
C
         C = DBLE(2*NVAR*NVAR+3*NVAR-1)/DBLE(6*(NVAR+1)*(NG-1))
         STAT = STAT*(1.0D0-C*(RSUM-1.0D0/WSUM))
         DF = 0.5D0*DBLE((NG-1)*NVAR*(NVAR+1))
         IFAULT = -1
         IF (STAT.GT.0.0D0) THEN
            SIG = G01ECF('U',STAT,DF,IFAULT)
         ELSE
            STAT = 0.0D0
            SIG = 1.0D0
         END IF
      END IF
  800 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, NVAR.lt.1: NVAR = ',I16)
99998 FORMAT (' ** On entry, NG.lt.2: NG = ',I16)
99997 FORMAT (' ** On entry, N.lt.1: N = ',I16)
99996 FORMAT (' ** On entry, M.lt.NVAR: M = ',I16,' NVAR = ',I16)
99995 FORMAT (' ** On entry, LDX.lt.N: LDX = ',I16,' N = ',I16)
99994 FORMAT (' ** On entry, LDG.lt.NG: LDG = ',I16,' NG = ',I16)
99993 FORMAT (' ** On entry, WEIGHT is not valid: WEIGHT = ',A1)
99992 FORMAT (' ** On entry, ',I16,' values of ISX.gt.0, not NVAR:',
     *       /'              NVAR = ',I16)
99991 FORMAT (' ** On entry, ',I16,' th value of ING not valid')
99990 FORMAT (' ** On entry, the ',I16,' th value of WT.lt.0.0')
99989 FORMAT (' ** The effective number of observations for group',I16,
     *       /'              is less than 1')
99988 FORMAT (' ** The number of observations for group',I16,/'       ',
     *       '      is less than NVAR')
99987 FORMAT (' ** The variables in group',I16,' are not of full rank')
99986 FORMAT (' ** The variables are not of full rank')
      END
