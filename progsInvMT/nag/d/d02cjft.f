      SUBROUTINE D02CJF(T,TEND,NEQ,Y,FCN,TOL,RELABS,OUTPUT,G,RWORK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED.  IER-680 (DEC 1989).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02CJF')
      DOUBLE PRECISION  RZERO, ONE, TEN, TENM7
      PARAMETER         (RZERO=0.0D0,ONE=1.0D0,TEN=10.0D0,TENM7=1.0D-7)
      INTEGER           MXREC
      PARAMETER         (MXREC=3)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TEND, TOL
      INTEGER           IFAIL, NEQ
      CHARACTER*1       RELABS
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(28+21*NEQ), Y(NEQ)
C     .. Function Arguments ..
      DOUBLE PRECISION  G
      EXTERNAL          G
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, OUTPUT
C     .. Scalars in Common ..
      DOUBLE PRECISION  DELSGN, FOURU, HOLD, SRBIG, SRU, SVTNEW, TLEFT,
     *                  TOLD, TROOTS, TSTAR, TWOU, U, U34, U78, XOLD,
     *                  XSAVE, ZERO
      INTEGER           IBEGIN, IINTEG, INFLOP, INIT, INROTP, IQUIT,
     *                  ITOL, ITSTOP, IVC, KGI, KLE4, KOLD, KORD, KPREV,
     *                  KROOTP, KSTEPS, NS
      LOGICAL           DISCOP, GSTOP, INTOUT, NEEDG, NEWGEQ, NORND,
     *                  PGSTOP, PHASE1, PSERCH, ROOT, ROOTS, SEARCH,
     *                  START, STIFF
      CHARACTER*6       GOPT, OUTOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  ALPHA(12), BETA(12), GI(11), GRWB(13), PSI(12),
     *                  SIG(13), V(12), W(12)
      INTEGER           IV(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  GDUM, GWANT, HMAX, TLAST, TOUT, TOUT1, TSAVE,
     *                  TSTOP, TWANT, UROUND
      INTEGER           BADCMP, I, IDID, IERR, IGNEW, IGOLD, IGP, IGSC,
     *                  IGV, IINDXG, IMMR, INEEDG, IP, IPHI, IPHI3N,
     *                  IPHI4N, IPHI5N, IPRO, IREVCM, IROD, ISEE, ITGV,
     *                  ITK, ITLB, ITRB, IWT, IYP, IYPOUT, IYY, KWANT,
     *                  MAXNUM, NEQG, NFAIL, NREC, NSUCC
      CHARACTER*1       RELAB1
C     .. Local Arrays ..
      DOUBLE PRECISION  ATOL(1), RTOL(1)
      INTEGER           IWORK(22)
      CHARACTER*80      REC(MXREC)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02QFX, E04UDU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, SQRT
C     .. Common blocks ..
      COMMON            /AD02CJ/OUTOPT, GOPT
      COMMON            /AD02QF/ALPHA, BETA, PSI, V, W, SIG, GRWB, GI,
     *                  XOLD, HOLD, TOLD, XSAVE, TSTAR, TWOU, INIT,
     *                  IBEGIN, ITOL, IINTEG, ITSTOP, INFLOP, IQUIT, IV,
     *                  NS, KORD, KOLD, KSTEPS, KLE4, KPREV, IVC, KGI,
     *                  START, PHASE1, NORND, STIFF, INTOUT
      COMMON            /BD02QF/ZERO, U, FOURU, SRU, U34, U78, SRBIG,
     *                  DELSGN, TROOTS, TLEFT, SVTNEW, KROOTP, INROTP,
     *                  GSTOP, PGSTOP, ROOT, ROOTS, NEEDG, DISCOP,
     *                  NEWGEQ, SEARCH, PSERCH
C     .. Save statement ..
      SAVE              /AD02CJ/, /AD02QF/, /BD02QF/
C     .. Executable Statements ..
      IERR = 0
      NREC = 0
      IF (NEQ.LT.1) THEN
         IERR = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99999) NEQ
      END IF
      IF (TOL.LE.RZERO) THEN
         IERR = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99998) TOL
      END IF
      IF (T.EQ.TEND) THEN
         IERR = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99997) T
      END IF
      IF (IERR.EQ.1) GO TO 220
      IDID = 0
      IBEGIN = 0
      ITOL = 0
      IINTEG = 0
      ITSTOP = 1
      TSTOP = TEND
      RWORK(1) = TEND
      IF (T.EQ.TEND) GO TO 220
      OUTOPT = 'YESOPT'
      TOUT = T
      CALL OUTPUT(TOUT,Y)
      IF (OUTOPT.EQ.'YESOPT') THEN
         IF (TOUT.EQ.T) THEN
            IERR = 4
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99988)
         ELSE IF (SIGN(ONE,TOUT-T).NE.SIGN(ONE,TEND-T)) THEN
            IERR = 4
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99987) TOUT
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99986) T, TEND
         END IF
         IF (IERR.EQ.4) GO TO 220
      ELSE
         TOUT = TEND
      END IF
      RTOL(1) = TOL
      ATOL(1) = TOL
      RELAB1 = RELABS
      CALL E04UDU(RELAB1)
      IF (RELAB1.EQ.'A') THEN
         RTOL(1) = RZERO
      ELSE IF (RELAB1.EQ.'R') THEN
         UROUND = X02AJF()
         ATOL(1) = MIN(SQRT(UROUND),MAX(TENM7,TEN*UROUND))
      ELSE IF (RELAB1.NE.'M' .AND. RELAB1.NE.'D') THEN
         IERR = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99996) RELAB1
         GO TO 220
      END IF
C
C     COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK ARRAY
C
C                             -- REAL RWORK
C
      IYPOUT = 21
      IYP = NEQ + IYPOUT
      IYY = NEQ + IYP
      IWT = NEQ + IYY
      IP = NEQ + IWT
      IPHI = NEQ + IP
C
      IPHI3N = IPHI + 3*NEQ
      IPHI4N = IPHI + 4*NEQ
      IPHI5N = IPHI + 5*NEQ
C
      GOPT = 'YSGOPT'
      GDUM = G(T,Y)
      IF (GOPT.EQ.'YSGOPT') THEN
         NEQG = 1
         ISEE = 0
         NEWGEQ = .TRUE.
         SEARCH = .FALSE.
         GSTOP = .TRUE.
      ELSE
         GSTOP = .FALSE.
         SEARCH = .FALSE.
         NEQG = 0
         ISEE = 0
      END IF
      IGOLD = 16*NEQ + IPHI
      IGNEW = NEQG + IGOLD
      IGP = NEQG + IGNEW
      ITK = NEQG + IGP
      IPRO = NEQG + ITK
      IROD = NEQG*ISEE + IPRO
      ITGV = NEQG*ISEE + IROD
      IGV = 3*NEQG*ISEE + ITGV
      ITLB = 3*NEQG*ISEE + IGV
      ITRB = NEQG*ISEE + ITLB
C
C
C                             -- INTEGER IWORK
      IGSC = 21
      IINDXG = NEQG*ISEE + IGSC
      IMMR = NEQG*ISEE + IINDXG
      INEEDG = NEQG*ISEE + IMMR
      IQUIT = 0
C
      TLAST = T
      HMAX = RZERO
      MAXNUM = 1000
      NSUCC = 0
      NFAIL = 0
      TSAVE = TOUT
      IF (OUTOPT.EQ.'YESOUT' .AND. ABS(TOUT-T).GT.ABS(TEND-T))
     *    TOUT = TEND
      IF (NEQG.EQ.0) NEQG = 1
C
      IREVCM = 0
   20 CONTINUE
      CALL D02QFX(IREVCM,TWANT,KWANT,GWANT,NEQ,T,Y,TOUT,RTOL,ATOL,IDID,
     *            RWORK(IYPOUT),RWORK(IYP),RWORK(IYY),RWORK(IWT),
     *            RWORK(IP),RWORK(IPHI),RWORK(IGOLD),RWORK(IGNEW),
     *            RWORK(ITGV),RWORK(IGV),RWORK(IGP),RWORK(ITK),
     *            RWORK(ITLB),RWORK(ITRB),RWORK(IPRO),RWORK(IROD),
     *            RWORK(1),RWORK(11),RWORK(12),RWORK(13),HMAX,MAXNUM,
     *            NSUCC,NFAIL,IWORK(IINDXG),IWORK(IGSC),IWORK(IMMR),
     *            IWORK(INEEDG),IWORK(11),IWORK(12),NEQG,BADCMP)
C
      GO TO (40,40,60,80,80,100,120,20,140,140,
     *       180,180) IREVCM
      GO TO 200
C
   40 CONTINUE
      CALL FCN(TWANT,Y,RWORK(IYPOUT))
      GO TO 20
   60 CONTINUE
      CALL FCN(TWANT,RWORK(IYY),RWORK(IPHI5N))
      GO TO 20
   80 CONTINUE
      CALL FCN(TWANT,RWORK(IPHI3N),RWORK(IPHI4N))
      GO TO 20
  100 CONTINUE
      CALL FCN(TWANT,RWORK(IP),RWORK(IYP))
      GO TO 20
  120 CONTINUE
      CALL FCN(TWANT,RWORK(IYY),RWORK(IYP))
      GO TO 20
  140 CONTINUE
      DO 160 I = 1, NEQG
         RWORK(IGOLD+I-1) = G(TWANT,Y)
  160 CONTINUE
      GO TO 20
  180 CONTINUE
      GWANT = G(TWANT,Y)
      GO TO 20
C
  200 CONTINUE
      IF (IDID.EQ.2 .OR. IDID.EQ.3) THEN
         IF (TOUT.NE.TEND) THEN
            TOUT1 = TOUT
            CALL OUTPUT(TOUT,Y)
            IF (TOUT.EQ.TOUT1) THEN
               IERR = 5
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99985) TOUT
            ELSE IF (SIGN(ONE,TOUT-TOUT1).NE.SIGN(ONE,TEND-TOUT1)) THEN
               IERR = 5
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99984) TOUT
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99983) TOUT1, TEND
            ELSE
               TSAVE = TOUT
               IF (ABS(TOUT-T).GT.ABS(TEND-T)) TOUT = TEND
               GO TO 20
            END IF
         ELSE IF (TSAVE.EQ.TEND) THEN
            CALL OUTPUT(TOUT,Y)
         END IF
         IF (GOPT.EQ.'YSGOPT') THEN
            IERR = 6
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99995)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99994)
         END IF
      ELSE IF (IDID.EQ.-1 .OR. IDID.EQ.-4) THEN
         IBEGIN = 1
         GO TO 20
      ELSE IF (IDID.EQ.-2) THEN
         IF (T.EQ.TLAST) THEN
            IERR = 3
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99993)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99992) TOL
         ELSE
            IERR = 2
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99991) T
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99990) TOL
         END IF
      ELSE IF (IDID.LT.3) THEN
         IERR = 7
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99989) IDID
      END IF
  220 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** N .lt. 1. N = ',I16,'.')
99998 FORMAT (' ** TOL .le. 0.0. TOL = ',1P,D13.5,'.')
99997 FORMAT (' ** X .eq. XEND. X = ',1P,D13.5,'.')
99996 FORMAT (' ** RELABS .NE. ''M'',''A'',''R'' or ''D''. RELABS = ',
     *       A1,'.')
99995 FORMAT (' ** No change in sign of the function G(T,Y) was detect',
     *       'ed ')
99994 FORMAT ('    in the integration range.')
99993 FORMAT (' ** No integration steps have been taken. Progress not ',
     *       'possible with the ')
99992 FORMAT ('    input value of TOL,',1P,D13.5,'.')
99991 FORMAT (' ** Integration successful as far as T = ',1P,D13.5,', ',
     *       'but further progress')
99990 FORMAT ('    not possible with the input value of TOL,',1P,D13.5,
     *       '.')
99989 FORMAT (' ** Impossible error - internal variable IDID = ',I16,
     *       '.')
99988 FORMAT (' ** On the first call to OUTPUT, XSOL remained unchange',
     *       'd.')
99987 FORMAT (' ** On the first call to OUTPUT, XSOL was returned as ',
     *       1P,D13.5,',')
99986 FORMAT ('    which is inconsistent with (X,XEND) - (',1P,D13.5,
     *       ',',D13.5,').')
99985 FORMAT (' ** On a call to OUTPUT, XSOL remained unchanged. XSOL ',
     *       '= ',1P,D13.5,'.')
99984 FORMAT (' ** On a call to OUTPUT, XSOL was returned as ',1P,D13.5,
     *       ', which is')
99983 FORMAT ('    inconsistent with previous XSOL and XEND - (',1P,
     *       D13.5,',',D13.5,').')
      END
