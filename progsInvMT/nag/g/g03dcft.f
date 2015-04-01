      SUBROUTINE G03DCF(TYPE,EQUAL,PRIORS,NVAR,NG,NIG,GMEAN,LDG,GC,DET,
     *                  NOBS,M,ISX,X,LDX,PRIOR,P,LDP,IAG,ATIQ,ATI,WK,
     *                  IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C      Allocates observations to groups according to posterior
C      probabilities.
C
C      for use after G03DAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03DCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDG, LDP, LDX, M, NG, NOBS, NVAR
      LOGICAL           ATIQ
      CHARACTER         EQUAL, PRIORS, TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  ATI(LDP,*), DET(NG), GC((NG+1)*NVAR*(NVAR+1)/2),
     *                  GMEAN(LDG,NVAR), P(LDP,NG), PRIOR(NG),
     *                  WK(2*NVAR), X(LDX,M)
      INTEGER           IAG(NOBS), ISX(M), NIG(NG)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, BA, BB, BP, BPDF, BQ, C, C1, CF, E, PIJ,
     *                  PMAX, RNIGT, SUM, TOL, UFLOW
      INTEGER           I, IERROR, IFAULT, IGP, IND, J, K, NREC, NU
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S14ABF, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          S14ABF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01EEF, G03DBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, DBLE
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 1
      IF (NVAR.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) NVAR
      ELSE IF (NG.LT.2) THEN
         WRITE (P01REC(1),FMT=99998) NG
      ELSE IF (NOBS.LT.1) THEN
         WRITE (P01REC(1),FMT=99997) NOBS
      ELSE IF (M.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99996) M, NVAR
      ELSE IF (LDG.LT.NG) THEN
         WRITE (P01REC(1),FMT=99995) LDG, NG
      ELSE IF (LDX.LT.NOBS) THEN
         WRITE (P01REC(1),FMT=99994) LDX, NOBS
      ELSE IF (LDP.LT.NOBS) THEN
         WRITE (P01REC(1),FMT=99993) LDP, NOBS
      ELSE IF (TYPE.NE.'E' .AND. TYPE.NE.'e' .AND. TYPE.NE.'P' .AND.
     *         TYPE.NE.'p') THEN
         WRITE (P01REC(1),FMT=99986) TYPE
      ELSE IF (EQUAL.NE.'E' .AND. EQUAL.NE.'e' .AND. EQUAL.NE.'U' .AND.
     *         EQUAL.NE.'u') THEN
         WRITE (P01REC(1),FMT=99992) EQUAL
      ELSE IF (PRIORS.NE.'E' .AND. PRIORS.NE.'e' .AND. PRIORS.NE.
     *         'I' .AND. PRIORS.NE.'i' .AND. PRIORS.NE.'P' .AND.
     *         PRIORS.NE.'p') THEN
         WRITE (P01REC(1),FMT=99991) PRIORS
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.1) GO TO 220
      UFLOW = LOG(X02AMF())
C
C     Check ISX
C
      K = 0
      DO 20 I = 1, M
         IF (ISX(I).GT.0) K = K + 1
   20 CONTINUE
      IF (K.NE.NVAR) THEN
         IERROR = 2
         WRITE (P01REC(1),FMT=99987) K, NVAR
         GO TO 220
      END IF
      NU = 0
      IF (EQUAL.EQ.'U' .OR. EQUAL.EQ.'u') THEN
         DO 40 I = 1, NG
            IF (NIG(I).LE.NVAR) THEN
               IERROR = 2
               WRITE (P01REC(1),FMT=99985) I
               GO TO 220
            END IF
            NU = NU + NIG(I)
   40    CONTINUE
      ELSE
         DO 60 I = 1, NG
            IF (NIG(I).LE.0) THEN
               IERROR = 2
               WRITE (P01REC(1),FMT=99983) I
               GO TO 220
            END IF
            NU = NU + NIG(I)
   60    CONTINUE
         IF (NU.LE.NG+NVAR) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99984)
            GO TO 220
         END IF
      END IF
      NU = NU - NG
      IND = 0
      IF (PRIORS.EQ.'I' .OR. PRIORS.EQ.'i') THEN
         IND = 1
         SUM = 0.0D0
         DO 80 I = 1, NG
            IF (PRIOR(I).LE.0.0D0) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99990) I
               GO TO 220
            END IF
            SUM = SUM + PRIOR(I)
   80    CONTINUE
         IF (ABS(SUM-1.0D0).GT.10.0D0*X02AJF()) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99989)
            GO TO 220
         END IF
      ELSE IF (PRIORS.EQ.'P' .OR. PRIORS.EQ.'p') THEN
         IND = 1
         RNIGT = DBLE(NU+NG)
         DO 100 I = 1, NG
            PRIOR(I) = DBLE(NIG(I))/RNIGT
  100    CONTINUE
      END IF
      IFAULT = 1
      CALL G03DBF(EQUAL,'S',NVAR,NG,GMEAN,LDG,GC,NOBS,M,ISX,X,LDX,P,LDP,
     *            WK,IFAULT)
      IF (IFAULT.NE.0) THEN
         IERROR = 4
         WRITE (P01REC(1),FMT=99988)
         GO TO 220
      END IF
      IF (EQUAL.EQ.'U' .OR. EQUAL.EQ.'u') IND = IND + 2
      IF (TYPE.EQ.'P' .OR. TYPE.EQ.'p') IND = IND + 4
C
C     IND now indicates which method is to be used
C
C                   equal cov                 unequal cov
C          equal priors  unequal priors equal priors  unequal priors
C
C     estim.    0               1            2               3
C     predict   4               5            6               7
C
      DO 140 J = 1, NG
C
C        Compute factors for estimative and predictive methods
C
         IFAULT = 0
         IF (IND.EQ.0 .OR. IND.EQ.4) THEN
            C1 = DBLE(NU+1)
            CF = DBLE(NIG(J))/DBLE(NU*(NIG(J)+1))
            A = 0.0D0
         ELSE IF (IND.EQ.1 .OR. IND.EQ.5) THEN
            C1 = DBLE(NU+1)
            CF = DBLE(NIG(J))/DBLE(NU*(NIG(J)+1))
            A = -2.0D0*LOG(PRIOR(J))
         ELSE IF (IND.EQ.2 .OR. IND.EQ.6) THEN
            C1 = DBLE(NIG(J))
            CF = C1/(C1*C1-1.0D0)
            A = DET(J)
         ELSE IF (IND.EQ.3 .OR. IND.EQ.7) THEN
            C1 = DBLE(NIG(J))
            CF = C1/(C1*C1-1.0D0)
            A = -2.0D0*LOG(PRIOR(J)) + DET(J)
         END IF
         IF (IND.GT.3) A = A - DBLE(NVAR)*LOG(CF)
         IF (IND.GT.5) A = A + 2.0D0*(S14ABF(0.5D0*(C1-DBLE(NVAR))
     *                     ,IFAULT)-S14ABF(0.5D0*C1,IFAULT))
C
C        set up constants for atypicality index
C
         IF (ATIQ) THEN
            TOL = 0.00005D0
            BA = 0.5D0*DBLE(NVAR)
            B = 1.0D0/CF
            BB = 0.5D0*C1 - BA
         END IF
C
C        Compute Posterior Probabilities (and atypicality index)
C
         DO 120 I = 1, NOBS
            PIJ = P(I,J)
            IF (ATIQ) THEN
               C = PIJ/(PIJ+B)
               IFAULT = -1
               CALL G01EEF(C,BA,BB,TOL,BP,BQ,BPDF,IFAULT)
               ATI(I,J) = BP
            END IF
            IF (IND.LE.3) THEN
               E = -0.5D0*(PIJ+A)
            ELSE
               E = -0.5D0*(C1*(LOG(1.0D0+CF*PIJ))+A)
            END IF
            IF (E.GT.UFLOW) THEN
               P(I,J) = EXP(E)
            ELSE
               P(I,J) = 0.0D0
            END IF
  120    CONTINUE
  140 CONTINUE
C
C     Standardize Posterior Probabilities and allocate to groups
C
      DO 200 I = 1, NOBS
         IGP = -1
         PMAX = 0.0D0
         SUM = 0.0D0
         DO 160 J = 1, NG
            IF (P(I,J).GT.PMAX) THEN
               PMAX = P(I,J)
               IGP = J
            END IF
            SUM = SUM + P(I,J)
  160    CONTINUE
         IAG(I) = IGP
         IF (SUM.GT.0.0D0) THEN
            DO 180 J = 1, NG
               P(I,J) = P(I,J)/SUM
  180       CONTINUE
         END IF
  200 CONTINUE
  220 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, NVAR.lt.1 : NVAR = ',I16)
99998 FORMAT (' ** On entry, NG.lt.2 : NG = ',I16)
99997 FORMAT (' ** On entry, NOBS.lt.1 : NOBS = ',I16)
99996 FORMAT (' ** On entry, M.lt.NVAR : M = ',I16,' NVAR = ',I16)
99995 FORMAT (' ** On entry, LDG.lt.NG : LDG = ',I16,' NG = ',I16)
99994 FORMAT (' ** On entry, LDX.lt.NOBS : LDX = ',I16,' NOBS = ',I16)
99993 FORMAT (' ** On entry, LDP.lt.NOBS : LDP = ',I16,' NOBS = ',I16)
99992 FORMAT (' ** On entry, EQUAL is not a valid character: EQUAL = ',
     *       A1)
99991 FORMAT (' ** On entry, PRIORS is not a valid character: PRIORS = '
     *       ,A1)
99990 FORMAT (' ** On entry, the ',I16,' th element of PRIOR.lt.0')
99989 FORMAT (' ** On entry, the sum of PRIOR.ne.1')
99988 FORMAT (' ** On entry, a diagonal element of R.eq.0')
99987 FORMAT (' ** On entry, ',I16,' values of ISX.gt.0, not NVAR = ',
     *       I16)
99986 FORMAT (' ** On entry, TYPE is not a valid character: TYPE = ',A1)
99985 FORMAT (' ** On entry, ',I16,' th value of NIG.le.NVAR')
99984 FORMAT (' ** On entry, the sum of NIG.le.NG+NVAR')
99983 FORMAT (' ** On entry, ',I16,' th value of NIG.le.0')
      END
