      SUBROUTINE G03ACF(WEIGHT,N,M,X,LDX,ISX,NX,ING,NG,WT,NIG,CVM,LDCVM,
     *                  E,LDE,NCV,CVX,LDCVX,TOL,IRANKX,WK,IWK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1671 (JUN 1995).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03ACF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, IRANKX, IWK, LDCVM, LDCVX, LDE, LDX, M,
     *                  N, NCV, NG, NX
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  CVM(LDCVM,NX), CVX(LDCVX,NG-1), E(LDE,6),
     *                  WK(IWK), WT(*), X(LDX,M)
      INTEGER           ING(N), ISX(M), NIG(NG)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, COND, E2, EPS, R, RDF, RN1, RN2, SCALE,
     *                  WSUM, WTOL
      INTEGER           I, IERROR, IFAULT, IRN1, IXPG, J, K, KK, MINDIM,
     *                  MINWK, MINXG, N0, NDV, NN1, NREC, NWK, NWT, NXN
      LOGICAL           SVDX, WTL
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP1(1,1), WKSP2(1,1), WKSP3(1,1)
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ, X02AJF
      INTEGER           F06KLF, P01ABF
      EXTERNAL          F02WDZ, X02AJF, F06KLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, DTRSV, F01QCF, F01QDF,
     *                  F02WEF, F02WUF, F06DBF, F06FBF, F06FCF, G03AAZ,
     *                  G03ACZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IXPG = NX + NG
      IF (NX.GE.NG-1) THEN
         NWK = (NX+1)*NX + (NX-1)*5
         MINXG = NG - 1
      ELSE
         NWK = (NX-1)*5 + (NG-1)*NX
         MINXG = NX
      END IF
      NN1 = N*NX
      NXN = N*NX + NX
      MINWK = MAX(N,NWK) + NN1
      IF (NX.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) NX
      ELSE IF (NG.LT.2) THEN
         WRITE (P01REC(1),FMT=99998) NG
      ELSE IF (M.LT.NX) THEN
         WRITE (P01REC(1),FMT=99997) M, NX
      ELSE IF (N.LT.IXPG) THEN
         NREC = 2
         WRITE (P01REC,FMT=99996) N, IXPG
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99995) LDX, N
      ELSE IF (LDCVX.LT.NX) THEN
         WRITE (P01REC(1),FMT=99994) LDCVX, NX
      ELSE IF (LDCVM.LT.NG) THEN
         WRITE (P01REC(1),FMT=99987) LDCVM, NG
      ELSE IF (LDE.LT.MINXG) THEN
         NREC = 2
         WRITE (P01REC,FMT=99993) LDE, MINXG
      ELSE IF (IWK.LT.MINWK) THEN
         NREC = 2
         WRITE (P01REC,FMT=99992) MINWK, IWK
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'V' .AND. WEIGHT.NE.'v' .AND. WEIGHT.NE.'U' .AND.
     *         WEIGHT.NE.'u') THEN
         WRITE (P01REC(1),FMT=99990) WEIGHT
      ELSE IF (TOL.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99984) TOL
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         EPS = X02AJF()
         IF (TOL.LT.EPS) THEN
            WTOL = SQRT(EPS)
         ELSE
            WTOL = TOL
         END IF
         CALL F06DBF(NG,0,NIG,1)
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w' .OR. WEIGHT.EQ.'V' .OR.
     *       WEIGHT.EQ.'v') THEN
            WTL = .TRUE.
            N0 = 0
            WSUM = 0.0D0
            NWT = (NX-1)*N
            DO 20 I = 1, N
               IF (WT(I).LT.0.0D0) THEN
                  IERROR = 2
                  WRITE (P01REC(1),FMT=99986) I
                  GO TO 460
               ELSE IF (WT(I).GT.0.0D0) THEN
                  IF (ING(I).LE.0 .OR. ING(I).GT.NG) GO TO 440
                  NIG(ING(I)) = NIG(ING(I)) + 1
                  N0 = N0 + 1
                  WSUM = WSUM + WT(I)
                  WK(NWT+I) = SQRT(WT(I))
               ELSE
                  WK(NWT+I) = 0.0D0
               END IF
   20       CONTINUE
         ELSE
            WTL = .FALSE.
            DO 40 I = 1, N
               IF (ING(I).LE.0 .OR. ING(I).GT.NG) GO TO 440
               NIG(ING(I)) = NIG(ING(I)) + 1
   40       CONTINUE
            N0 = N
            WSUM = DBLE(N)
         END IF
         K = 0
         DO 60 I = 1, M
            IF (ISX(I).GT.0) K = K + 1
   60    CONTINUE
         IF (K.NE.NX) THEN
            IERROR = 4
            WRITE (P01REC(1),FMT=99989) K, NX
         ELSE
            IF (WTL) THEN
               CALL G03AAZ('U','W',N,X,LDX,M,ISX,NX,WT,WSUM,WK,N,WKSP1,
     *                     WK(NN1+1))
            ELSE
               CALL G03AAZ('U','U',N,X,LDX,M,ISX,NX,WT,WSUM,WK,N,WKSP1,
     *                     WK(NN1+1))
            END IF
C
C           CALCULATE GROUP MEANS FOR QX
C
            CALL F01QCF(N,NX,WK,N,CVX,IFAULT)
            K = 0
            DO 120 I = 1, NG
               IF (NIG(I).NE.0) THEN
                  K = K + 1
                  IF (WTL) THEN
                     DO 80 J = 1, N
                        IF (ING(J).EQ.I) THEN
                           WK(NN1+J) = SQRT(WT(J))
                        ELSE
                           WK(NN1+J) = 0.0D0
                        END IF
   80                CONTINUE
                  ELSE
                     DO 100 J = 1, N
                        IF (ING(J).EQ.I) THEN
                           WK(NN1+J) = 1.0D0
                        ELSE
                           WK(NN1+J) = 0.0D0
                        END IF
  100                CONTINUE
                  END IF
                  IFAULT = 1
                  CALL F01QDF('T','S',N,NX,WK,N,CVX,1,WK(NN1+1),N,WKSP1,
     *                        IFAULT)
                  CALL DCOPY(NX,WK(NN1+1),1,CVM(I,1),LDCVM)
               ELSE
                  CALL F06FBF(NX,0.0D0,CVM(I,1),LDCVM)
               END IF
  120       CONTINUE
            NDV = K - 1
            IF (NDV.LE.0) THEN
               IERROR = 7
               WRITE (P01REC(1),FMT=99983)
               GO TO 460
            ELSE IF (N0-(NDV+NX).LE.0) THEN
               IERROR = 7
               NREC = 2
               WRITE (P01REC,FMT=99981)
               GO TO 460
            END IF
            IF (WEIGHT.EQ.'V' .OR. WEIGHT.EQ.'v') WSUM = DBLE(N0)
            RDF = SQRT(WSUM-DBLE(K))
            COND = F02WDZ(NX,WK,N,WK(NN1+1))
            IF (COND*WTOL.GT.1.0D0) THEN
               SVDX = .TRUE.
               IFAULT = 1
               CALL F02WUF(NX,WK,N,0,WKSP1,1,.TRUE.,WK(NXN+1),NX,
     *                     WK(NN1+1),.TRUE.,WK(NXN+NX*NX+1),IFAULT)
               IF (IFAULT.GT.0) THEN
                  IERROR = 5
                  WRITE (P01REC(1),FMT=99991)
                  GO TO 460
C
               ELSE
                  IRANKX = F06KLF(NX,WK(NN1+1),1,WTOL)
                  IF (IRANKX.LE.0) THEN
                     IERROR = 8
                     WRITE (P01REC(1),FMT=99982)
                     GO TO 460
                  END IF
                  DO 140 I = 1, IRANKX
                     WK(NN1+I) = 1.0D0/WK(NN1+I)
  140             CONTINUE
                  DO 160 I = 1, NX
                     CALL F06FCF(IRANKX,WK(NN1+1),1,WK((I-1)*N+1),1)
  160             CONTINUE
                  DO 180 I = 1, NG
                     IF (NIG(I).NE.0) THEN
                        CALL DCOPY(NX,CVM(I,1),LDCVM,WK(NN1+1),1)
                        CALL DGEMV('T',NX,IRANKX,1.0D0,WK(NXN+1),NX,
     *                             WK(NN1+1),1,0.0D0,CVM(I,1),LDCVM)
                     END IF
  180             CONTINUE
               END IF
            ELSE
               SVDX = .FALSE.
               IRANKX = NX
            END IF
            K = 0
  200       CONTINUE
            K = K + 1
            IRN1 = NIG(K)
            IF (IRN1.EQ.0) GO TO 200
            IF (WTL) THEN
               CALL F06FBF(NG-1,0.0D0,CVX(IRANKX,1),LDCVX)
               RN1 = 0.0D0
               DO 220 I = 1, N
                  IF (WT(I).GT.0.0D0) THEN
                     KK = ING(I)
                     IF (KK.EQ.K) THEN
                        RN1 = RN1 + WT(I)
                     ELSE
                        KK = KK - 1
                        CVX(IRANKX,KK) = CVX(IRANKX,KK) + WT(I)
                     END IF
                  END IF
  220          CONTINUE
            ELSE
               RN1 = DBLE(IRN1)
               DO 240 I = K + 1, NG
                  CVX(IRANKX,I-1) = DBLE(NIG(I))
  240          CONTINUE
            END IF
            CALL DCOPY(IRANKX,CVM(K,1),LDCVM,WK(NN1+1),1)
            SCALE = 1.0D0/RN1
            CALL DSCAL(IRANKX,SCALE,CVM(K,1),LDCVM)
            KK = 0
            DO 280 I = K + 1, NG
               IF (NIG(I).NE.0) THEN
                  KK = KK + 1
                  RN2 = CVX(IRANKX,I-1)
                  SCALE = 1.0D0/RN2
                  R = RN1*SCALE
                  A = 1.0D0/SQRT((1.0D0+R)*RN1)
                  B = -R*A
                  DO 260 J = 1, IRANKX
                     CVX(J,KK) = WK(NN1+J)*A + CVM(I,J)*B
                     WK(NN1+J) = WK(NN1+J) + CVM(I,J)
                     CVM(I,J) = CVM(I,J)*SCALE
  260             CONTINUE
                  RN1 = RN1 + RN2
               END IF
  280       CONTINUE
            IFAULT = 1
            IF (IRANKX.GE.NDV) THEN
               MINDIM = NDV
               CALL F02WEF(IRANKX,NDV,CVX,LDCVX,0,WKSP1,1,.TRUE.,WKSP2,
     *                     1,E(1,1),.FALSE.,WKSP3,1,WK(NN1+1),IFAULT)
               IF (IFAULT.GT.0) THEN
                  IERROR = 5
                  WRITE (P01REC(1),FMT=99991)
                  GO TO 460
C
               END IF
            ELSE
               MINDIM = IRANKX
               DO 300 I = 1, NDV
                  CALL DCOPY(IRANKX,CVX(1,I),1,WK(NN1+(I-1)*IRANKX+1),1)
  300          CONTINUE
               CALL F02WEF(IRANKX,NDV,WK(NN1+1),IRANKX,0,WKSP1,1,.TRUE.,
     *                     CVX,LDCVX,E(1,1),.FALSE.,WKSP2,1,
     *                     WK(NN1+NDV*NX+1),IFAULT)
               IF (IFAULT.GT.0) THEN
                  IERROR = 5
                  WRITE (P01REC(1),FMT=99991)
                  GO TO 460
C
               END IF
            END IF
            NCV = F06KLF(MINDIM,E(1,1),1,WTOL)
            DO 320 I = 1, NCV
               E2 = 1.0D0 - E(I,1)*E(I,1)
               IF (E2.LE.0.0D0) GO TO 340
               SCALE = RDF/SQRT(E2)
               CALL DSCAL(IRANKX,SCALE,CVX(1,I),1)
  320       CONTINUE
            GO TO 360
  340       IERROR = 6
            WRITE (P01REC,FMT=99985)
            NCV = 1
  360       CONTINUE
            DO 380 I = 1, NG
               IF (NIG(I).NE.0) THEN
                  CALL DCOPY(IRANKX,CVM(I,1),LDCVM,WK(NN1+1),1)
                  CALL DGEMV('T',IRANKX,NCV,1.0D0,CVX,LDCVX,WK(NN1+1),1,
     *                       0.0D0,CVM(I,1),LDCVM)
               END IF
  380       CONTINUE
            IF (SVDX) THEN
               DO 400 I = 1, NCV
                  CALL DCOPY(IRANKX,CVX(1,I),1,WK(NN1+1),1)
                  CALL DGEMV('T',IRANKX,NX,1.0D0,WK,N,WK(NN1+1),1,0.0D0,
     *                       CVX(1,I),1)
  400          CONTINUE
            ELSE
               DO 420 I = 1, NCV
                  CALL DTRSV('U','N','N',NX,WK,N,CVX(1,I),1)
  420          CONTINUE
            END IF
C
C           CALCULATE TEST STATISTICS
C
            CALL G03ACZ(E,LDE,WSUM,NCV,IRANKX,NDV,IFAULT)
         END IF
         GO TO 460
C
  440    IERROR = 3
         WRITE (P01REC(1),FMT=99988) I
      END IF
  460 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry NX .lt. 1: NX = ',I16)
99998 FORMAT (' ** On entry NG .lt. 2: NG = ',I16)
99997 FORMAT (' ** On entry M .lt. NX: M = ',I16,' NX = ',I16)
99996 FORMAT (' ** On entry N .lt. NX+NG: N = ',I16,/'                ',
     *       '       NX+NG = ',I16)
99995 FORMAT (' ** On entry LDX .lt. N: LDX = ',I16,' N = ',I16)
99994 FORMAT (' ** On entry LDCVX .lt. NX: LDCVX = ',I16,' NX = ',I16)
99993 FORMAT (' ** On entry LDE .lt. MIN(NX,NG-1): LDE = ',I16,/'     ',
     *       '                      MIN(NX,NG-1) = ',I16)
99992 FORMAT (' ** On entry IWK is too small, min value = ',I16,/'    ',
     *       '                                 IWK = ',I16)
99991 FORMAT (' ** An SVD has failed to converge')
99990 FORMAT (' ** On entry WEIGHT is not valid: WEIGHT = ',A1)
99989 FORMAT (' ** On entry there are ',I16,' X vars instead of ',I16)
99988 FORMAT (' ** On entry ',I16,' th value of ING not valid')
99987 FORMAT (' ** On entry LDCVM .lt. NG: LDCVM = ',I16,' NG = ',I16)
99986 FORMAT (' ** On entry the ',I16,' th value of WT .lt. 0.0')
99985 FORMAT (' ** Canonical correlation equal to 1.0')
99984 FORMAT (' ** On entry TOL .lt. 0.0: TOL = ',D13.5)
99983 FORMAT (' ** Less than 2 groups have non-zero membership')
99982 FORMAT (' ** The rank of X is 0')
99981 FORMAT (' ** The effective number of observations is less than',/
     *     '    the effective number of groups plus number of variables'
     *       )
      END
