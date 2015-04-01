      SUBROUTINE G02EAF(MEAN,WEIGHT,N,M,X,LDX,NAME,ISX,Y,WT,NMOD,MODEL,
     *                  LDM,RSS,NTERMS,MRANK,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES ALL POSSIBLE REGRESION USING QR METHOD
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDM, LDX, M, N, NMOD
      CHARACTER         MEAN, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  RSS(LDM), WK(N*(M+1)), WT(*), X(LDX,M), Y(N)
      INTEGER           ISX(M), MRANK(LDM), NTERMS(LDM)
      CHARACTER*(*)     MODEL(LDM,M)
      CHARACTER*(*)     NAME(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, COND, RSSF, S, SUM, TEMP, TOL, WTSUM
      INTEGER           I, IERROR, IFAULT, IFR, IK, ISS, ISWAP, IT, J,
     *                  JJ, K, KF, L, LK, MS, MT, N0, NN, NN1, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ, X02AJF
      INTEGER           G02EAZ, P01ABF
      EXTERNAL          F02WDZ, X02AJF, G02EAZ, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01QCF, F06BAF, M01DAF, M01DBF, M01EAF, M01ECF,
     *                  DCOPY, DROT, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (N.LT.2) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (M.LT.2) THEN
         WRITE (P01REC(1),FMT=99998) M
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99997) LDX, N
      ELSE IF (LDM.LT.M) THEN
         WRITE (P01REC(1),FMT=99989) LDM, M
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         IF ((MEAN.NE.'Z' .AND. MEAN.NE.'M')
     *       .AND. (MEAN.NE.'z' .AND. MEAN.NE.'m')) THEN
            IERROR = 1
            WRITE (P01REC(1),FMT=99996) MEAN
         ELSE IF ((WEIGHT.NE.'U' .AND. WEIGHT.NE.'W')
     *            .AND. (WEIGHT.NE.'u' .AND. WEIGHT.NE.'w')) THEN
            IERROR = 1
            WRITE (P01REC(1),FMT=99995) WEIGHT
         ELSE
            TOL = 10.0D0*X02AJF()
            IK = 0
C
C           SELECT FREE AND KERNAL VARIABLES
C
            DO 20 J = 1, M
               IF (ISX(J).LT.0) THEN
                  IERROR = 3
                  WRITE (P01REC(1),FMT=99988) J
                  GO TO 800
               ELSE IF (ISX(J).GE.2) THEN
                  IK = IK + 1
                  MRANK(IK) = J
               END IF
   20       CONTINUE
            IT = IK
            DO 40 J = 1, M
               IF (ISX(J).EQ.1) THEN
                  IT = IT + 1
                  MRANK(IT) = J
               END IF
   40       CONTINUE
            IFR = IT - IK
            IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
               MT = IT + 1
            ELSE
               MT = IT
            END IF
            IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
               IF (MT.GE.N) THEN
                  IERROR = 5
                  WRITE (P01REC,FMT=99991)
                  GO TO 800
               END IF
            END IF
            IF (IFR.EQ.0) THEN
               IERROR = 3
               WRITE (P01REC(1),FMT=99992)
            ELSE
               NMOD = 2**IFR
               IF (LDM.LT.NMOD) THEN
                  IERROR = 4
                  NREC = 2
                  WRITE (P01REC,FMT=99990) LDM, NMOD
               ELSE
C
C                 STANDARDIZE AND TRANSFORM X VARIABLES IF REQUIRED
C
                  IF ((WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u')
     *                .AND. (MEAN.EQ.'z' .OR. MEAN.EQ.'Z')) THEN
                     N0 = N
                     DO 60 J = 1, IT
                        CALL DCOPY(N,X(1,MRANK(J)),1,WK((J-1)*N+1),1)
   60                CONTINUE
                  ELSE IF ((WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u')
     *                     .AND. (MEAN.EQ.'M' .OR. MEAN.EQ.'m')) THEN
                     N0 = N
                     DO 120 J = 1, IT
                        JJ = MRANK(J)
                        SUM = 0.0D0
                        DO 80 I = 1, N
                           SUM = X(I,JJ) + SUM
   80                   CONTINUE
                        SUM = SUM/N
                        DO 100 I = 1, N
                           WK((J-1)*N+I) = X(I,JJ) - SUM
  100                   CONTINUE
  120                CONTINUE
                  ELSE IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
C
C                    CHECK WEIGHTS
C
                     N0 = 0
                     WTSUM = 0.0D0
                     DO 140 I = 1, N
                        IF (WT(I).LT.0.0D0) THEN
                           GO TO 260
C
                        ELSE IF (WT(I).GT.0.0D0) THEN
                           N0 = N0 + 1
                           WTSUM = WTSUM + WT(I)
                           WK(IT*N+I) = SQRT(WT(I))
                        ELSE
                           WK(IT*N+I) = 0.0D0
                        END IF
  140                CONTINUE
                     IF (MT.GE.N0) THEN
                        IERROR = 5
                        WRITE (P01REC,FMT=99991)
                        GO TO 800
                     END IF
                     IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
                        DO 180 J = 1, IT
                           DO 160 I = 1, N
                              WK((J-1)*N+I) = X(I,MRANK(J))*WK(IT*N+I)
  160                      CONTINUE
  180                   CONTINUE
                     ELSE
                        DO 240 J = 1, IT
                           JJ = MRANK(J)
                           SUM = 0.0D0
                           DO 200 I = 1, N
                              SUM = WT(I)*X(I,JJ) + SUM
  200                      CONTINUE
                           SUM = SUM/WTSUM
                           DO 220 I = 1, N
                              WK((J-1)*N+I) = WK(IT*N+I)*(X(I,MRANK(J))
     *                                        -SUM)
  220                      CONTINUE
  240                   CONTINUE
                     END IF
                     GO TO 280
C
  260                IERROR = 2
                     WRITE (P01REC(1),FMT=99994) I
                     GO TO 800
C
                  END IF
  280             CONTINUE
C
C                 STANDARDIZE AND TRANSFORM Y
C
                  IF ((WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u')
     *                .AND. (MEAN.EQ.'z' .OR. MEAN.EQ.'Z')) THEN
                     CALL DCOPY(N,Y,1,WK(N*IT+1),1)
                  ELSE IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
                     SUM = 0.0D0
                     DO 300 I = 1, N
                        SUM = Y(I) + SUM
  300                CONTINUE
                     SUM = SUM/N0
                     DO 320 I = 1, N
                        WK(IT*N+I) = Y(I) - SUM
  320                CONTINUE
                  ELSE IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
                     DO 340 I = 1, N
                        WK(IT*N+I) = WK(IT*N+I)*Y(I)
  340                CONTINUE
                  ELSE
                     SUM = 0.0D0
                     DO 360 I = 1, N
                        SUM = WT(I)*Y(I) + SUM
  360                CONTINUE
                     SUM = SUM/WTSUM
                     DO 380 I = 1, N
                        WK(IT*N+I) = WK(IT*N+I)*(Y(I)-SUM)
  380                CONTINUE
                  END IF
C
C                 PERFORM QR
C
                  IFAULT = 1
                  CALL F01QCF(N,IT+1,WK,N,RSS,IFAULT)
                  COND = F02WDZ(IT,WK,N,RSS)
                  IF (COND*TOL.GT.1.0D0) THEN
                     IERROR = 6
                     WRITE (P01REC,FMT=99993)
                  ELSE
                     MS = IFR*IFR
C
C                    CONTRACT STORAGE SPACE FOR R
C
                     DO 420 I = 1, IFR
                        DO 400 J = 1, I
                           WK((I-1)*IFR+J) = WK((IK+I-1)*N+IK+J)
  400                   CONTINUE
  420                CONTINUE
                     DO 440 I = 1, IFR
                        WK(MS+I) = WK(IT*N+IK+I)
  440                CONTINUE
                     RSSF = WK(N*IT+IT+1)
                     DO 480 J = 1, IK
                        K = MRANK(J)
                        DO 460 I = 1, NMOD
                           MODEL(I,J) = NAME(K)
  460                   CONTINUE
  480                CONTINUE
C
C                    CALCULATE RSS FROM INITIAL FULL MODEL
C
                     RSSF = RSSF*RSSF
                     SUM = 0.0D0
                     NN = NMOD - IFR
                     ISS = NN + 1
                     RSS(NN) = RSSF
                     NTERMS(NN) = IFR + IK
                     DO 500 I = 1, IFR
                        MODEL(NN,I+IK) = NAME(MRANK(I+IK))
  500                CONTINUE
                     DO 540 I = IFR, 1, -1
                        NN = NN + 1
                        SUM = WK(MS+I)**2 + SUM
                        RSS(NN) = RSSF + SUM
                        DO 520 L = 1, I - 1
                           MODEL(NN,L+IK) = NAME(MRANK(L+IK))
  520                   CONTINUE
  540                CONTINUE
C
C                    GENERATE MODEL SEQUENCE
C
                     IF (IFR.GT.1) THEN
                        L = 1
                        NTERMS(1) = 1
                        DO 600 I = 3, IFR
                           LK = L
                           DO 560 J = I - 1, 1, -1
                              L = L + 1
                              NTERMS(L) = J
  560                      CONTINUE
                           DO 580 J = 1, LK
                              L = L + 1
                              NTERMS(L) = NTERMS(J) + 1
  580                      CONTINUE
  600                   CONTINUE
C
C                       CALCULATE RSS FOR MODEL SEQUENCE
C
                        NN = NMOD - IFR - 1
                        DO 680 KF = IFR, 1, -1
                           NN1 = 2**(KF-1) - KF + 1
                           DO 660 I = NN, NN1, -1
                              J = NTERMS(I)
                              ISWAP = MRANK(J+IK)
                              MRANK(J+IK) = MRANK(J+IK+1)
                              MRANK(J+IK+1) = ISWAP
                              CALL M01DBF(MRANK(IK+1),1,J,'A',
     *                                    NTERMS(ISS),IFAULT)
                              DO 620 K = 1, J
                                 MODEL(I,NTERMS(ISS-1+K)+IK)
     *                             = NAME(MRANK(K+IK))
  620                         CONTINUE
                              CALL F06BAF(WK(IFR*J+J),WK(IFR*J+J+1),C,S)
                              WK(IFR*J+J+1) = 0.0D0
                              CALL DSWAP(J,WK(IFR*(J-1)+1),1,WK(IFR*J+1)
     *                                   ,1)
                              CALL DROT(KF-J,WK(IFR*J+J),IFR,
     *                                  WK(IFR*J+J+1),IFR,C,S)
                              TEMP = WK(MS+J)
                              WK(MS+J) = WK(MS+J+1)*S + C*TEMP
                              WK(MS+J+1) = WK(MS+J+1)*C - S*TEMP
                              SUM = 0.0D0
                              DO 640 L = J + 1, IFR
                                 SUM = WK(MS+L)**2 + SUM
  640                         CONTINUE
                              RSS(I) = SUM + RSSF
                              NTERMS(I) = J + IK
  660                      CONTINUE
                           NN = NN1 - 1
  680                   CONTINUE
                     END IF
                     K = IT
                     DO 700 I = ISS, NMOD
                        K = K - 1
                        NTERMS(I) = K
  700                CONTINUE
                     IFAULT = 1
                     CALL M01DBF(NTERMS(1),1,NMOD,'A',MRANK,IFAULT)
                     IFAULT = 1
                     CALL M01EAF(RSS,1,NMOD,MRANK,IFAULT)
                     DO 720 J = 1, M
                        IFAULT = 1
                        CALL M01ECF(MODEL(1,J),1,NMOD,MRANK,IFAULT)
  720                CONTINUE
                     NTERMS(1) = IK
                     NN = 2
                     NN1 = 1
                     DO 780 K = 1, IFR - 1
                        IFAULT = 1
                        NN1 = G02EAZ(K,IFR,IFAULT) + NN1
                        DO 740 I = NN, NN1
                           NTERMS(I) = K + IK
  740                   CONTINUE
                        IFAULT = 1
                        CALL M01DAF(RSS,NN,NN1,'D',MRANK,IFAULT)
                        IFAULT = 1
                        CALL M01EAF(RSS,NN,NN1,MRANK,IFAULT)
                        DO 760 J = 1, M
                           IFAULT = 1
                           CALL M01ECF(MODEL(1,J),NN,NN1,MRANK,IFAULT)
  760                   CONTINUE
                        NN = NN1 + 1
  780                CONTINUE
                     NTERMS(NMOD) = IT
                     IFAULT = 1
                     CALL M01DAF(RSS,1,NMOD,'A',MRANK,IFAULT)
                  END IF
               END IF
            END IF
         END IF
      END IF
  800 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, N.lt.2 : N = ',I16)
99998 FORMAT (' ** On entry, M.lt.2 : M = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N : LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, MEAN is not valid : MEAN = ',A1)
99995 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99994 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
99993 FORMAT (' ** Full model is not of full rank')
99992 FORMAT (' ** There are no free X variables')
99991 FORMAT (' ** Number of requested x-variables .ge. number',' of o',
     *       'bservations')
99990 FORMAT (' ** On entry, LDM is too small : LDM = ',I16,/'    it m',
     *       'ust be at least = ',I16)
99989 FORMAT (' ** On entry, LDM.lt.M : LDM = ',I16,' M = ',I16)
99988 FORMAT (' ** On entry, ISX(',I16,').lt.0')
      END
