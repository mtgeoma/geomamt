      SUBROUTINE G02DAF(MEAN,WEIGHT,N,X,LDX,M,ISX,IP,Y,WT,RSS,IDF,B,SE,
     *                  COV,RES,H,Q,LDQ,SVD,IRANK,P,TOL,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15B REVISED. IER-954 (NOV 1991).
C     MARK 17 REVISED. IER-1658 (JUN 1995).
C     COMPUTES LINEAR REGRESSION
C
C     The non-full rank case is covered.
C     A QR/SVD approach is used.
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02DAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS, TOL
      INTEGER           IDF, IFAIL, IP, IRANK, LDQ, LDX, M, N
      LOGICAL           SVD
      CHARACTER         MEAN, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV(IP*(IP+1)/2), H(N), P(IP*(IP+2)),
     *                  Q(LDQ,IP+1), RES(N), SE(IP), WK(5*(IP-1)+IP*IP),
     *                  WT(N), X(LDX,M), Y(N)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, RMS, SQWT
      INTEGER           I, IERROR, IFAULT, IJ, IP2, J, K, NO, NREC
      LOGICAL           MEANL, WTL
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP(1)
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ, DDOT
      INTEGER           F06KLF, P01ABF
      EXTERNAL          F02WDZ, DDOT, F06KLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01QCF, F01QDF, F02WUF, F06FBF, F06FCF, G02AAX,
     *                  G02AAY, G02AAZ, DCOPY, DGEMV, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      SVD = .FALSE.
      NREC = 1
      IERROR = 1
      IF (N.LT.2) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (M.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) M
      ELSE IF (IP.LT.1) THEN
         WRITE (P01REC(1),FMT=99988) IP
      ELSE IF (IP.GT.N) THEN
         WRITE (P01REC(1),FMT=99992) IP, N
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99997) LDX, N
      ELSE IF (LDQ.LT.N) THEN
         WRITE (P01REC(1),FMT=99989) LDQ, N
      ELSE IF (TOL.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99990) TOL
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
            MEANL = .FALSE.
         ELSE IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            MEANL = .TRUE.
         ELSE
            IERROR = 2
            WRITE (P01REC(1),FMT=99996) MEAN
            GO TO 360
C
         END IF
         IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
            WTL = .FALSE.
         ELSE IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            WTL = .TRUE.
         ELSE
            IERROR = 2
            WRITE (P01REC(1),FMT=99995) WEIGHT
            GO TO 360
C
         END IF
C
C        CHECK WEIGHTS
C
         IF (WTL) THEN
            NO = 0
            DO 20 I = 1, N
               IF (WT(I).LT.0.0D0) THEN
                  GO TO 40
C
               ELSE IF (WT(I).GT.0.0D0) THEN
                  NO = NO + 1
                  RES(I) = SQRT(WT(I))
               ELSE
                  RES(I) = 0.0D0
               END IF
   20       CONTINUE
            GO TO 60
C
   40       IERROR = 3
            WRITE (P01REC(1),FMT=99994) I
            GO TO 360
C
         ELSE
C
C
            NO = N
         END IF
   60    CALL F06FBF(IP*(IP+1)/2,0.0D0,COV,1)
         IF (MEANL) THEN
            K = 1
         ELSE
            K = 0
         END IF
         DO 80 J = 1, M
            IF (ISX(J).LT.0) THEN
               IERROR = 4
               WRITE (P01REC(1),FMT=99986) J
               GO TO 360
            ELSE IF (ISX(J).GT.0) THEN
               K = K + 1
            END IF
   80    CONTINUE
         IF (IP.NE.K) THEN
            IERROR = 4
            NREC = 2
            WRITE (P01REC,FMT=99987) IP
         ELSE
            IF (WTL) THEN
               IF (MEANL) THEN
                  K = 1
                  CALL DCOPY(N,RES,1,Q(1,2),1)
               ELSE
                  K = 0
               END IF
               DO 120 J = 1, M
                  IF (ISX(J).LT.0) THEN
                     IERROR = 4
                     WRITE (P01REC(1),FMT=99986) J
                     GO TO 360
                  ELSE IF (ISX(J).GT.0) THEN
                     K = K + 1
                     DO 100 I = 1, N
                        Q(I,K+1) = X(I,J)*RES(I)
  100                CONTINUE
                  END IF
  120          CONTINUE
               DO 140 I = 1, N
                  Q(I,1) = RES(I)*Y(I)
  140          CONTINUE
            ELSE
               IF (MEANL) THEN
                  K = 1
                  CALL F06FBF(N,1.0D0,Q(1,2),1)
               ELSE
                  K = 0
               END IF
               DO 160 J = 1, M
                  IF (ISX(J).LT.0) THEN
                     IERROR = 4
                     WRITE (P01REC(1),FMT=99986) J
                     GO TO 360
                  ELSE IF (ISX(J).GT.0) THEN
                     K = K + 1
                     CALL DCOPY(N,X(1,J),1,Q(1,K+1),1)
                  END IF
  160          CONTINUE
               CALL DCOPY(N,Y,1,Q,1)
            END IF
            IFAULT = 1
            CALL F01QCF(N,IP,Q(1,2),LDQ,P,IFAULT)
            COND = F02WDZ(IP,Q(1,2),LDQ,WK)
            IFAULT = 1
            CALL F01QDF('T','S',N,IP,Q(1,2),LDQ,P,1,Q,N,WKSP,IFAULT)
            IF (COND*TOL.GT.1.0D0) THEN
               SVD = .TRUE.
               DO 180 I = 2, IP + 1
                  CALL DCOPY(I-1,Q(1,I),1,P(I*IP+1),1)
  180          CONTINUE
               CALL DCOPY(IP,Q,1,SE,1)
               IP2 = IP*IP
               IFAULT = 1
               CALL F02WUF(IP,P(2*IP+1),IP,1,SE,IP,.TRUE.,WK,IP,P(IP+1),
     *                     .TRUE.,WK(IP2+1),IFAULT)
               IF (IFAULT.LT.0) THEN
                  IERROR = 6
                  WRITE (P01REC(1),FMT=99993)
               ELSE
                  IRANK = F06KLF(IP,P(IP+1),1,TOL)
                  DO 200 I = 1, IRANK
                     WK(IP2+I) = 1.0D0/P(IP+I)
  200             CONTINUE
                  DO 220 I = 2, IP + 1
                     CALL F06FCF(IRANK,WK(IP2+1),1,P(I*IP+1),1)
  220             CONTINUE
                  CALL DGEMV('T',IRANK,IP,1.0D0,P(2*IP+1),IP,SE,1,0.0D0,
     *                       B,1)
                  IJ = 1
                  DO 240 I = 1, IP
                     CALL DCOPY(IRANK,P(IP*I+IP+1),1,WK(IP2+1),1)
                     CALL DGEMV('T',IRANK,I,1.0D0,P(2*IP+1),IP,WK(IP2+1)
     *                          ,1,0.0D0,COV(IJ),1)
                     IJ = IJ + I
  240             CONTINUE
                  DO 260 I = 1, N
                     IF(WTL) THEN
                        SQWT = SQRT(WT(I))
                        IF (MEANL) THEN
                           RES(1) = SQWT
                           K = 1
                        ELSE
                           K = 0
                        END IF
                        DO 245 J = 1, M
                           IF (ISX(J).GT.0) THEN
                              K = K + 1
                              RES(K) = SQWT*X(I,J)
                           END IF
 245                    CONTINUE
                     ELSE
                        IF (MEANL) THEN
                           RES(1) = 1.0D0
                           K = 1
                        ELSE
                           K = 0
                        END IF
                        DO 250 J = 1, M
                           IF (ISX(J).GT.0) THEN
                              K = K + 1
                              RES(K) = X(I,J)
                           END IF
  250                   CONTINUE
                     END IF
                     CALL DGEMV('N',IRANK,IP,1.0D0,P(2*IP+1),IP,RES,
     *                             1,0.0D0,WK(IP2+1),1)
                     H(I) = DDOT(IRANK,WK(IP2+1),1,WK(IP2+1),1)
  260             CONTINUE
                  IDF = NO - IRANK
                  IF (IDF.EQ.0) THEN
                     IERROR = 5
                     WRITE (P01REC(1),FMT=99991)
                     RSS = 0
                     CALL F06FBF(N,0.0D0,RES,1)
                     GO TO 360
                  END IF
                  IF (IP.GT.IRANK) THEN
                     CALL DGEMV('N',IP,IP-IRANK,1.0D0,WK(IRANK*IP+1),IP,
     *                          SE(IRANK+1),1,0.0D0,RES,1)
                  ELSE
                     CALL F06FBF(IP,0.0D0,RES,1)
                  END IF
                  CALL DCOPY(N-IP,Q(IP+1,1),1,RES(IP+1),1)
                  IFAULT = 1
                  CALL F01QDF('N','S',N,IP,Q(1,2),LDQ,P,1,RES,N,WKSP,
     *                        IFAULT)
               END IF
            ELSE
               IRANK = IP
               CALL DCOPY(IP,Q,1,B,1)
               CALL DTRSV('U','N','N',IP,Q(1,2),LDQ,B,1)
               IJ = 1
               DO 280 I = 1, IP
                  CALL DCOPY(I,Q(1,I+1),1,COV(IJ),1)
                  IJ = IJ + I
  280          CONTINUE
               CALL G02AAZ('U','N',IP,COV)
               CALL G02AAX('U',IP,COV,WK)
               CALL G02AAY('L','N',IP,WK)
               CALL G02AAX('L',IP,WK,COV)
               IFAULT = 1
               DO 300 I = 1, N
                   IF(WTL) THEN
                     SQWT = SQRT(WT(I))
                     IF (MEANL) THEN
                        RES(1) = SQWT
                        K = 1
                     ELSE
                        K = 0
                     END IF
                     DO 305 J = 1, M
                        IF (ISX(J).GT.0) THEN
                           K = K + 1
                           RES(K) = SQWT*X(I,J)            
                        END IF
 305                 CONTINUE
                   ELSE
                       IF (MEANL) THEN
                        RES(1) = 1.0D0
                        K = 1
                     ELSE
                        K = 0
                     END IF
                     DO 310 J = 1, M
                        IF (ISX(J).GT.0) THEN
                           K = K + 1
                           RES(K) = X(I,J)            
                         END IF
 310                 CONTINUE
                   END IF
                   CALL DTRSV('U','T','N',IP,Q(1,2),LDQ,RES,1)
                   H(I) = DDOT(IP,RES,1,RES,1)
  300          CONTINUE
               IDF = NO - IRANK
               IF (IDF.EQ.0) THEN
                  IERROR = 5
                  WRITE (P01REC(1),FMT=99991)
                  CALL F06FBF(N,0.0D0,RES,1)
                  RSS = 0
                  GO TO 360
               END IF
               CALL F06FBF(IP,0.0D0,RES,1)
               CALL DCOPY(N-IP,Q(IP+1,1),1,RES(IP+1),1)
               IFAULT = 1
               CALL F01QDF('N','S',N,IP,Q(1,2),LDQ,P,1,RES,N,WKSP,
     *                     IFAULT)
            END IF
            RSS = DDOT(N,RES,1,RES,1)
            RMS = RSS/DBLE(IDF)
            DO 320 I = 1, (IP*IP+IP)/2
               COV(I) = COV(I)*RMS
  320       CONTINUE
            DO 340 J = 1, IP
               SE(J) = SQRT(COV((J*J+J)/2))
  340       CONTINUE
         END IF
      END IF
  360 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry N .lt. 2 : N = ',I16)
99998 FORMAT (' ** On entry M .lt. 1 : M = ',I16)
99997 FORMAT (' ** On entry LDX .lt. N : LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry MEAN is not a valid character : MEAN is ',
     *       A1)
99995 FORMAT (' ** On entry WEIGHT is not a valid character :WEIGHT is '
     *       ,A1)
99994 FORMAT (' ** On entry the ',I16,' th value of WT .lt. 0')
99993 FORMAT (' ** SVD solution failed to converge')
99992 FORMAT (' ** On entry IP .gt. N : IP = ',I16,' N = ',I16)
99991 FORMAT (' ** Degrees of freedom for error equal 0')
99990 FORMAT (' ** On entry TOL .lt. 0.0: TOL = ',D12.5)
99989 FORMAT (' ** On entry LDQ .lt. N : LDQ = ',I16,' N = ',I16)
99988 FORMAT (' ** On entry IP .lt. 1 : IP = ',I16)
99987 FORMAT (' ** On entry value of IP incompatible with number of no',
     *       'n-zero values of ISX:',/'             IP = ',I16)
99986 FORMAT (' ** On entry the ',I16,' th value of ISX .lt. 0')
      END
