      SUBROUTINE G02EEF(ISTEP,MEAN,WEIGHT,N,M,X,LDX,NAME,ISX,MAXIP,Y,WT,
     *                  FIN,ADDVAR,NEWVAR,CHRSS,F,MODEL,NTERM,RSS,IDF,
     *                  IFR,FREE,EXSS,Q,LDQ,P,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     FITS A REGRESSION MODEL BY FORWARD SELECTION
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02EEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHRSS, F, FIN, RSS
      INTEGER           IDF, IFAIL, IFR, ISTEP, LDQ, LDX, M, MAXIP, N,
     *                  NTERM
      LOGICAL           ADDVAR
      CHARACTER*1       MEAN, WEIGHT
      CHARACTER*(*)     NEWVAR
C     .. Array Arguments ..
      DOUBLE PRECISION  EXSS(MAXIP), P(MAXIP+1), Q(LDQ,MAXIP+2),
     *                  WK(2*MAXIP), WT(*), X(LDX,M), Y(N)
      INTEGER           ISX(M)
      CHARACTER*(*)     FREE(MAXIP)
      CHARACTER*(*)     MODEL(MAXIP)
      CHARACTER*(*)     NAME(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, COND, QY, RSSN, TOL, XQY, XSS, ZETA
      INTEGER           I, IERROR, IFAULT, IMEAN, IN, IP, J, MAXP, NO,
     *                  NREC
      LOGICAL           MEANL, WTL
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ, DDOT, X02AJF
      INTEGER           P01ABF
      EXTERNAL          F02WDZ, DDOT, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01QCF, F01QDF, F06FBF, F06FCF, F06FTF, G02AAT,
     *                  G02EEZ, DCOPY, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      ADDVAR = .FALSE.
      IF (N.LT.2) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (M.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) M
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99997) LDX, N
      ELSE IF (LDQ.LT.N) THEN
         WRITE (P01REC(1),FMT=99988) LDQ, N
      ELSE IF (ISTEP.LT.0) THEN
         WRITE (P01REC(1),FMT=99993) ISTEP
      ELSE IF (ISTEP.NE.0 .AND. NTERM.LE.0) THEN
         NREC = 2
         WRITE (P01REC,FMT=99992) ISTEP, NTERM
      ELSE IF (ISTEP.NE.0 .AND. RSS.LE.0.0D0) THEN
         WRITE (P01REC(1),FMT=99985) RSS
      ELSE IF (FIN.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99983) FIN
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         TOL = X02AJF()
         IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
            MEANL = .FALSE.
            IMEAN = 0
         ELSE IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            MEANL = .TRUE.
            IMEAN = 1
         ELSE
            IERROR = 1
            WRITE (P01REC(1),FMT=99996) MEAN
            GO TO 180
C
         END IF
         IF (ISTEP.EQ.0) THEN
C
C        SET UP Q MATRIX IF FIRST CALL
C        PICK KERNAL AND FREE VARIABLES
C
            IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
               WTL = .FALSE.
            ELSE IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
               WTL = .TRUE.
            ELSE
               IERROR = 1
               WRITE (P01REC(1),FMT=99995) WEIGHT
               GO TO 180
C
            END IF
C
C     CHECK WEIGHTS
C
            IF (WTL) THEN
               NO = 0
               DO 20 I = 1, N
                  IF (WT(I).LT.0.0D0) THEN
                     GO TO 40
C
                  ELSE IF (WT(I).GT.0.0D0) THEN
                     NO = NO + 1
                     Q(I,1) = SQRT(WT(I))
                  ELSE
                     Q(I,1) = 0.0D0
                  END IF
   20          CONTINUE
               GO TO 60
C
   40          IERROR = 2
               WRITE (P01REC(1),FMT=99994) I
               GO TO 180
C
            ELSE
               NO = N
            END IF
   60       NTERM = 0
            IF (MEANL .AND. WTL) THEN
               CALL DCOPY(N,Q,1,Q(1,2),1)
            ELSE IF (MEANL) THEN
               CALL F06FBF(N,1.0D0,Q(1,2),1)
            END IF
            DO 80 J = 1, M
               IF (ISX(J).LT.0) THEN
                  IERROR = 4
                  WRITE (P01REC(1),FMT=99982) J
                  GO TO 180
               ELSE IF (ISX(J).GE.2) THEN
                  NTERM = NTERM + 1
                  IF (NTERM.GT.MAXIP) THEN
                     IERROR = 4
                     NREC = 2
                     WRITE (P01REC,FMT=99981) MAXIP
                     GO TO 180
                  END IF
                  MODEL(NTERM) = NAME(J)
                  IF (WTL) THEN
                     CALL G02AAT(N,Q,1,X(1,J),1,Q(1,NTERM+IMEAN+1),1)
                  ELSE
                     CALL DCOPY(N,X(1,J),1,Q(1,NTERM+IMEAN+1),1)
                  END IF
               END IF
   80       CONTINUE
            MAXP = NTERM + IMEAN
            IDF = NO - MAXP
            IF (IDF.LE.0) THEN
               IERROR = 3
               WRITE (P01REC(1),FMT=99984)
               GO TO 180
C
            ELSE
               IFR = 0
               DO 100 J = 1, M
                  IF (ISX(J).EQ.1) THEN
                     MAXP = MAXP + 1
                     IF (MAXP-IMEAN.GT.MAXIP) THEN
                        IERROR = 4
                        NREC = 2
                        WRITE (P01REC,FMT=99981) MAXIP
                        GO TO 180
                     END IF
                     IFR = IFR + 1
                     FREE(IFR) = NAME(J)
                     IF (WTL) THEN
                        CALL G02AAT(N,Q,1,X(1,J),1,Q(1,MAXP+1),1)
                     ELSE
                        CALL DCOPY(N,X(1,J),1,Q(1,MAXP+1),1)
                     END IF
                  END IF
  100          CONTINUE
               IF (MAXP.EQ.0) THEN
                  IERROR = 4
                  WRITE (P01REC(1),FMT=99990)
                  GO TO 180
C
               ELSE
                  IF (WTL) THEN
                     CALL F06FCF(N,Y,1,Q,1)
                  ELSE
                     CALL DCOPY(N,Y,1,Q,1)
                  END IF
                  IP = NTERM + IMEAN
                  IF (IP.GE.1) THEN
C
C      ADD KERNAL TERMS TO MODEL
C
                     IFAULT = 1
                     CALL F01QCF(N,IP,Q(1,2),LDQ,P,IFAULT)
                     COND = F02WDZ(IP,Q(1,2),LDQ,WK)
                     IFAULT = 1
                     CALL F01QDF('T','S',N,IP,Q(1,2),LDQ,P,IFR,Q(1,IP+2)
     *                           ,LDQ,WK,IFAULT)
                     CALL F01QDF('T','S',N,IP,Q(1,2),LDQ,P,1,Q,N,WK,
     *                           IFAULT)
                     IF (COND*TOL.GT.1.0D0) THEN
                        IERROR = 5
                        WRITE (P01REC(1),FMT=99989)
                        GO TO 180
                     END IF
                  END IF
                  RSS = DDOT(N-IP,Q(IP+1,1),1,Q(IP+1,1),1)
               END IF
            END IF
         END IF
         ISTEP = ISTEP + 1
         IF (IDF.LT.1) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99991)
         ELSE IF (IFR.EQ.0) THEN
            IERROR = 6
            WRITE (P01REC(1),FMT=99987)
         ELSE
C
C     SELECT BEST TERM
C
            CHRSS = 0.0D0
            IP = NTERM + IMEAN + 1
            QY = Q(IP,1)
            DO 120 I = 1, IFR
               CALL G02EEZ(0,N-IP,Q(IP,IP+I),Q(IP+1,IP+I),1,TOL,WK(I),
     *                     ZETA)
               IF (ZETA.EQ.0) THEN
                  EXSS(I) = 0.0D0
               ELSE
                  XQY = DDOT(N-IP,Q(IP+1,1),1,Q(IP+1,IP+I),1)
                  XSS = (1.0D0-ZETA*ZETA)*QY + XQY/WK(I)
                  EXSS(I) = XSS*XSS
                  WK(IFR+I) = ZETA
                  IF (EXSS(I).GT.CHRSS) THEN
                     CHRSS = EXSS(I)
                     IN = I
                  END IF
               END IF
  120       CONTINUE
            RSSN = RSS - CHRSS
            IF (RSSN.LE.0.0D0) THEN
               IERROR = 7
               WRITE (P01REC,FMT=99986)
            ELSE
               F = (IDF-1)*CHRSS/RSSN
               IF (F.GT.FIN) THEN
C
C      ADD NEW TERM IF SIGNIFICANT
C
                  ADDVAR = .TRUE.
                  NTERM = NTERM + 1
                  NEWVAR = FREE(IN)
                  IF (IN.NE.1) THEN
                     FREE(IN) = FREE(1)
                     EXSS(IN) = EXSS(1)
                  END IF
                  BETA = WK(IN)
                  ZETA = WK(IFR+IN)
                  IN = IN + IP - 1
                  MODEL(NTERM) = NEWVAR
                  IFR = IFR - 1
                  IDF = IDF - 1
                  CALL G02EEZ(1,N-IP,Q(IP,IN+1),Q(IP+1,IN+1),1,TOL,BETA,
     *                        ZETA)
                  IF (IN.NE.IP) CALL DSWAP(N,Q(1,IP+1),1,Q(1,IN+1),1)
                  CALL F06FTF(N-IP,Q(IP,1),Q(IP+1,1),1,ZETA,Q(IP+1,IP+1)
     *                        ,1)
                  RSS = DDOT(N-IP,Q(IP+1,1),1,Q(IP+1,1),1)
                  DO 140 I = IP + 2, IP + IFR + 1
                     CALL F06FTF(N-IP,Q(IP,I),Q(IP+1,I),1,ZETA,
     *                           Q(IP+1,IP+1),1)
  140             CONTINUE
                  P(IP) = ZETA
                  DO 160 I = 1, IFR
                     FREE(I) = FREE(I+1)
                     EXSS(I) = EXSS(I+1)
  160             CONTINUE
               END IF
            END IF
         END IF
      END IF
  180 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, N.lt.2 : N = ',I16)
99998 FORMAT (' ** On entry, M.lt.1 : M = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N : LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, MEAN is not valid : MEAN = ',A1)
99995 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99994 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
99993 FORMAT (' ** On entry, ISTEP.lt.0 : ISTEP = ',I16)
99992 FORMAT (' ** On entry, ISTEP AND NTERM are inconsistent:',/'    ',
     *       '         ISTEP = ',I16,' NTERM = ',I16)
99991 FORMAT (' ** Degrees of freedom for error will equal 0 if new va',
     *       'riable is  added')
99990 FORMAT (' ** Maximum number of variables to be included is 0')
99989 FORMAT (' ** Forced variables not of full rank')
99988 FORMAT (' ** On entry, LDQ.lt.N : LDQ = ',I16,' N = ',I16)
99987 FORMAT (' ** There are no free variables in the regression')
99986 FORMAT (' ** Denominator of F statistic is .le. 0.0')
99985 FORMAT (' ** On entry, with non-zero ISTEP, RSS.le.0.0 : RSS = ',
     *       D13.5)
99984 FORMAT (' ** On entry, number of forced variables .ge. N, i.e. I',
     *       'DF would be zero.')
99983 FORMAT (' ** On entry, FIN.lt.0.0 : FIN = ',D13.5)
99982 FORMAT (' ** On entry, ISX(',I16,').lt.0')
99981 FORMAT (' ** On entry, MAXIP is too small for number of terms gi',
     *       'ven by ISX:',/'                  MAXIP = ',I16)
      END
