      SUBROUTINE D02QGF(NEQF,T,Y,TOUT,NEQG,ROOT,IREVCM,TRVCM,YRVCM,
     *                  YPRVCM,GRVCM,KGRVCM,RWORK,LRWORK,IWORK,LIWORK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 17 REVISED. IER-1632 (JUN 1995).
C     .. Parameters ..
      INTEGER           ASK, OKAY, SET
      PARAMETER         (ASK=0,OKAY=1,SET=1)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02QGF')
      INTEGER           SVNEQF, SVNEQG, SVALTG, SV1STP, SVSPHS, SVLRWK,
     *                  SVLIWK, SVCRIT, SVVCTL, SVSTAT, SVNSUC, SVNFAI,
     *                  SVROOT, SVBADC, SVMXST, SVKOLD, SVKORD, CRITSV,
     *                  HMAXSV, HOLDSV, HNXTSV, TLFCSV, TCURSV
      PARAMETER         (SVNEQF=1,SVNEQG=2,SVALTG=3,SV1STP=4,SVSPHS=5,
     *                  SVLRWK=6,SVLIWK=7,SVCRIT=8,SVVCTL=9,SVSTAT=10,
     *                  SVNSUC=13,SVNFAI=14,SVROOT=15,SVBADC=16,
     *                  SVMXST=17,SVKOLD=18,SVKORD=19,CRITSV=1,HMAXSV=2,
     *                  HOLDSV=10,HNXTSV=11,TLFCSV=12,TCURSV=13)
      INTEGER           MAXREC
      PARAMETER         (MAXREC=10)
      DOUBLE PRECISION  RZERO, ONE
      PARAMETER         (RZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GRVCM, T, TOUT, TRVCM
      INTEGER           IFAIL, IREVCM, KGRVCM, LIWORK, LRWORK, NEQF,
     *                  NEQG, YPRVCM, YRVCM
      LOGICAL           ROOT
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), Y(NEQF)
      INTEGER           IWORK(LIWORK)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DELSGN, FOURU, HOLD, SRBIG, SRU, SVTNEW, TCHECK,
     *                  TLEFT, TOLD, TROOTS, TSTAR, TWOU, U, U34, U78,
     *                  XOLD, XSAVE, ZERO
      INTEGER           IATOL, IBEGIN, IGNEW, IGOLD, IGP, IGSC, IGV,
     *                  IINDXG, IINTEG, IMMR, INEEDG, INFLOP, INIT,
     *                  INROTP, IP, IPHI, IPHI3N, IPHI4N, IPHI5N, IPRO,
     *                  IQUIT, IROD, ITGV, ITK, ITLB, ITOL, ITRB,
     *                  ITSTOP, IVC, IWT, IYP, IYPOUT, IYY, KGI, KLE4,
     *                  KOLD, KORD, KPREV, KROOTP, KSTEPS, NEQGCP, NS
      LOGICAL           DISCOP, GSTOP, INTOUT, LROOT, NEEDG, NEWGEQ,
     *                  NORND, PGSTOP, PHASE1, PSERCH, ROOTS, SEARCH,
     *                  START, STIFF
C     .. Arrays in Common ..
      DOUBLE PRECISION  ALPHA(12), BETA(12), GI(11), GRWB(13), PSI(12),
     *                  SIG(13), V(12), W(12)
      INTEGER           IV(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  TSTOP
      INTEGER           ERSTAT, IDID, IER, IRTOL, ISEE, NREC, STATE
      LOGICAL           LOOP, VECTOL
C     .. Local Arrays ..
      CHARACTER*80      REC(MAXREC)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02QFX, D02QFY, D02QWZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN
C     .. Common blocks ..
      COMMON            /AD02QF/ALPHA, BETA, PSI, V, W, SIG, GRWB, GI,
     *                  XOLD, HOLD, TOLD, XSAVE, TSTAR, TWOU, INIT,
     *                  IBEGIN, ITOL, IINTEG, ITSTOP, INFLOP, IQUIT, IV,
     *                  NS, KORD, KOLD, KSTEPS, KLE4, KPREV, IVC, KGI,
     *                  START, PHASE1, NORND, STIFF, INTOUT
      COMMON            /BD02QF/ZERO, U, FOURU, SRU, U34, U78, SRBIG,
     *                  DELSGN, TROOTS, TLEFT, SVTNEW, KROOTP, INROTP,
     *                  GSTOP, PGSTOP, LROOT, ROOTS, NEEDG, DISCOP,
     *                  NEWGEQ, SEARCH, PSERCH
      COMMON            /MD02QF/TCHECK
CRWB - added IRTOL to the following SAVE'd COMMON block since needed
      COMMON            /ND02QF/IATOL, IGNEW, IGOLD, IGP, IGSC, IGV,
     *                  IINDXG, IMMR, INEEDG, IP, IPHI, IPHI3N, IPHI4N,
     *                  IPHI5N, IPRO, IROD, ITGV, ITK, ITLB, ITRB, IWT,
     *                  IYP, IYPOUT, IYY, NEQGCP, IRTOL
C     .. Save statement ..
      SAVE              /AD02QF/, /BD02QF/, /MD02QF/, /ND02QF/, LOOP
C     .. Data statements ..
      DATA              LOOP/.FALSE./
C     .. Executable Statements ..
C
C     JUMP ON REVERSE COMMUNICATION CALL
C
      GO TO (20,20,20,20,20,20,20,20,160,
     *       160,20,20) IREVCM
C
      ROOT = .FALSE.
      IF (IREVCM.NE.0) THEN
         NREC = 1
         IER = 1
         WRITE (REC(1),FMT=99971) IREVCM
         GO TO 220
      END IF
C
C     CHECK PREVIOUS CALL TO SETUP ROUTINE
C
      CALL D02QWZ(STATE,ASK)
      IF (STATE.NE.OKAY) THEN
         IER = 1
         IF (STATE.EQ.0) THEN
            NREC = 1
            WRITE (REC(1),FMT=99999)
         ELSE
            NREC = 2
            WRITE (REC(1),FMT=99998)
            WRITE (REC(2),FMT=99997) STATE - 1
         END IF
         GO TO 220
      END IF
C
      IER = 0
      NREC = 0
C
C     SET CONTROL VARIABLES FORM SETUP
C
      IBEGIN = IWORK(SVSTAT)
      ITOL = IWORK(SVVCTL)
      IINTEG = IWORK(SV1STP)
      ITSTOP = IWORK(SVCRIT)
C
C     SET SOME OUTPUT QUANTITIES
C
      IWORK(SVROOT) = 0
      IWORK(SVBADC) = 0
      IF (IWORK(SVSTAT).EQ.0) THEN
         RWORK(HNXTSV) = RZERO
         HOLD = RZERO
         KOLD = 0
         KORD = 0
         RWORK(HOLDSV) = HOLD
         RWORK(TLFCSV) = ONE
         RWORK(TCURSV) = T
         IWORK(SVNSUC) = 0
         IWORK(SVNFAI) = 0
      END IF
C
C     DATA CHECKS
C
      IF (T.EQ.TOUT) THEN
         IER = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99984) T
      END IF
      IF (IWORK(SVLRWK).NE.LRWORK) THEN
         IER = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99996)
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99995) IWORK(SVLRWK), LRWORK
      END IF
      IF (IWORK(SVLIWK).NE.LIWORK) THEN
         IER = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99994)
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99993) IWORK(SVLIWK), LIWORK
      END IF
      IF (IWORK(SVNEQF).NE.NEQF) THEN
         IER = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99992)
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99991) IWORK(SVNEQF), NEQF
      END IF
      IF (NEQG.LE.0) THEN
         IF (IWORK(SVNEQG).GT.0) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99990) NEQG
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99989) IWORK(SVNEQG)
         END IF
      ELSE
         IF (IWORK(SVNEQG).EQ.0) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99988) NEQG
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99987)
         ELSE IF (NEQG.NE.IWORK(SVNEQG)) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99986) NEQG
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99985) IWORK(SVNEQG)
         END IF
      END IF
      IF (IER.EQ.1) GO TO 220
C
C     SET WORKSPACE POINTERS
C
      IYPOUT = 21
      IYP = NEQF + IYPOUT
      IYY = NEQF + IYP
      IWT = NEQF + IYY
      IP = NEQF + IWT
      IPHI = NEQF + IP
      IPHI3N = IPHI + 3*NEQF
      IPHI4N = IPHI + 4*NEQF
      IPHI5N = IPHI + 5*NEQF
      IRTOL = 16*NEQF + IPHI
      VECTOL = IWORK(SVVCTL) .EQ. 1
      IF (VECTOL) THEN
         IATOL = NEQF + IRTOL
      ELSE
         IATOL = 1 + IRTOL
      END IF
C
C     CHECK FOR ROOTFINDING CAPABILITY
C
      IF (IWORK(SVNEQG).EQ.0) THEN
         GSTOP = .FALSE.
         SEARCH = .FALSE.
         ISEE = 0
         NEQGCP = 0
      ELSE
         GSTOP = .TRUE.
         NEQGCP = NEQG
         NEWGEQ = IWORK(SVALTG) .EQ. 1
         SEARCH = IWORK(SVSPHS) .EQ. 1
         IF (SEARCH) THEN
            ISEE = 1
         ELSE
            ISEE = 0
         END IF
      END IF
C
C     SET REMAINING WORKSPACE POINTERS
C
      IF (VECTOL) THEN
         IGOLD = NEQF + IATOL
      ELSE
         IGOLD = 1 + IATOL
      END IF
      IGNEW = NEQGCP + IGOLD
      IGP = NEQGCP + IGNEW
      ITK = NEQGCP + IGP
      IPRO = NEQGCP + ITK
      IROD = NEQGCP*ISEE + IPRO
      ITGV = NEQGCP*ISEE + IROD
      IGV = 3*NEQGCP*ISEE + ITGV
      ITLB = 3*NEQGCP*ISEE + IGV
      ITRB = NEQGCP*ISEE + ITLB
C
      IGSC = 21
      IINDXG = NEQGCP*ISEE + IGSC
      IMMR = NEQGCP*ISEE + IINDXG
      INEEDG = NEQGCP*ISEE + IMMR
C
C     CHECK DATA ON CONTINUATION CALL
C
      IF (IBEGIN.EQ.1) THEN
         IF (T.NE.TOLD) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99981) TOLD
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99980) T
         ELSE IF (INIT.NE.1) THEN
            IF (DELSGN*(TOUT-T).LT.RZERO) THEN
               IER = 1
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99979) TOUT
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99978)
            END IF
         END IF
         IF (IER.EQ.1) GO TO 220
      END IF
C
C     CHECK STOPPING POINT
C
      IF (ITSTOP.EQ.1) THEN
         TSTOP = RWORK(CRITSV)
         IF ( .NOT. (SIGN(ONE,TOUT-T).EQ.SIGN(ONE,TSTOP-T)
     *       .AND. ABS(TOUT-T).LE.ABS(TSTOP-T))) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99983) TOUT
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99982) TSTOP
            GO TO 220
         END IF
      END IF
C
      TSTAR = T
      IF (NEQGCP.EQ.0) NEQGCP = 1
      IREVCM = 0
      IDID = 0
   20 CONTINUE
      CALL D02QFX(IREVCM,TRVCM,KGRVCM,GRVCM,NEQF,T,Y,TOUT,RWORK(IRTOL),
     *            RWORK(IATOL),IDID,RWORK(IYPOUT),RWORK(IYP),RWORK(IYY),
     *            RWORK(IWT),RWORK(IP),RWORK(IPHI),RWORK(IGOLD),
     *            RWORK(IGNEW),RWORK(ITGV),RWORK(IGV),RWORK(IGP),
     *            RWORK(ITK),RWORK(ITLB),RWORK(ITRB),RWORK(IPRO),
     *            RWORK(IROD),RWORK(CRITSV),RWORK(HNXTSV),RWORK(TLFCSV),
     *            RWORK(TCURSV),RWORK(HMAXSV),IWORK(SVMXST),
     *            IWORK(SVNSUC),IWORK(SVNFAI),IWORK(IINDXG),IWORK(IGSC),
     *            IWORK(IMMR),IWORK(INEEDG),IWORK(11),IWORK(12),NEQGCP,
     *            IWORK(SVBADC))
C
      GO TO (40,40,60,80,80,100,120,240,140,
     *       140,180,180) IREVCM
      GO TO 200
C
   40 CONTINUE
C      CALL F(NEQF,TRVCM,Y,RWORK(IYPOUT))
C      GO TO 500
      YRVCM = 0
      YPRVCM = IYPOUT
      GO TO 240
   60 CONTINUE
C      CALL F(NEQF,TRVCM,RWORK(IYY),RWORK(IPHI5N))
C      GO TO 500
      YRVCM = IYY
      YPRVCM = IPHI5N
      GO TO 240
   80 CONTINUE
C      CALL F(NEQF,TRVCM,RWORK(IPHI3N),RWORK(IPHI4N))
C      GO TO 500
      YRVCM = IPHI3N
      YPRVCM = IPHI4N
      GO TO 240
  100 CONTINUE
C      CALL F(NEQF,TRVCM,RWORK(IP),RWORK(IYP))
C      GO TO 500
      YRVCM = IP
      YPRVCM = IYP
      GO TO 240
  120 CONTINUE
C      CALL F(NEQF,TRVCM,RWORK(IYY),RWORK(IYP))
C      GO TO 500
      YRVCM = IYY
      YPRVCM = IYP
      GO TO 240
C
C
  140 CONTINUE
C      DO 640 I = 1, NEQGCP
C         RWORK(IGOLD+I-1) = G(NEQF,TRVCM,Y,RWORK(IYPOUT),I)
C     640 CONTINUE
      KGRVCM = 1
      YRVCM = 0
      YPRVCM = IYPOUT
      GO TO 240
  160 RWORK(IGOLD+KGRVCM-1) = GRVCM
      KGRVCM = KGRVCM + 1
      IF (KGRVCM.LE.NEQGCP) THEN
         GO TO 240
      ELSE
         GO TO 20
      END IF
  180 CONTINUE
C      GRVCM = G(NEQF,TRVCM,Y,RWORK(IYPOUT),KGRVCM)
      YRVCM = 0
      YPRVCM = IYPOUT
      GO TO 240
C
  200 CONTINUE
      XSAVE = RWORK(TCURSV)
      RWORK(HOLDSV) = HOLD
      IWORK(SVSTAT) = IBEGIN
      IWORK(SVKOLD) = KOLD
      IWORK(SVKORD) = KORD
      IF (IDID.GT.0) THEN
         IER = 0
         ROOT = IDID .GT. 3
         IF (ROOT) IWORK(SVROOT) = 1
      ELSE IF (IDID.EQ.-1) THEN
         IER = 2
         NREC = 1
         WRITE (REC(1),FMT=99977)
      ELSE IF (IDID.EQ.-2) THEN
         IER = 3
         NREC = 1
         WRITE (REC(1),FMT=99976)
      ELSE IF (IDID.EQ.-3) THEN
         IER = 4
         NREC = 1
         WRITE (REC(1),FMT=99975)
      ELSE IF (IDID.EQ.-4) THEN
         IER = 5
         NREC = 1
         WRITE (REC(1),FMT=99974)
      ELSE IF (IDID.EQ.-8) THEN
         IER = 6
         NREC = 1
         WRITE (REC(1),FMT=99973)
      END IF
C
  220 CONTINUE
C
      IF (IER.GT.0) THEN
         IF ( .NOT. LOOP) THEN
            LOOP = .TRUE.
            TCHECK = T
         ELSE
            IF (TCHECK.EQ.T) THEN
               IER = 7
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99972) T
            ELSE
               LOOP = .FALSE.
            END IF
         END IF
      ELSE
         LOOP = .FALSE.
      END IF
C
      IREVCM = 0
      ERSTAT = IER + 1
      CALL D02QFY(ERSTAT,SET)
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
  240 IFAIL = 0
      RETURN
C
99999 FORMAT (' ** The setup routine D02QWF has not been called.')
99998 FORMAT (' ** The previous call to the setup routine D02QWF resul',
     *       'ted in the ')
99997 FORMAT ('    error exit IFAIL = ',I2,'.')
99996 FORMAT (' ** The dimension of RWORK supplied to D02QGF is not th',
     *       'at supplied to D02QWF:')
99995 FORMAT ('    LRWORK in D02QWF was ',I16,', LRWORK in D02QGF is ',
     *       I16,'.')
99994 FORMAT (' ** The dimension of IWORK supplied to D02QGF is not th',
     *       'at supplied to D02QWF:')
99993 FORMAT ('    LIWORK in D02QWF was ',I16,', LIWORK in D02QGF is ',
     *       I16,'.')
99992 FORMAT (' ** The value of NEQF supplied to D02QGF is not that su',
     *       'pplied to D02QWF:')
99991 FORMAT ('    NEQF in D02QWF was ',I16,', NEQF in D02QGF is ',I16,
     *       '.')
99990 FORMAT (' ** The non-positive value of NEQG,',I16,' indicating r',
     *       'ootfinding not required')
99989 FORMAT ('    conflicts with the positive value of NEQG,',I16,', ',
     *       'supplied to D02QWF.')
99988 FORMAT (' ** The positive value of NEQG,',I16,', indicating root',
     *       'finding required')
99987 FORMAT ('    conflicts with the non-positive value of NEQG suppl',
     *       'ied to D02QWF.')
99986 FORMAT (' ** The positive value of NEQG,',I16,', indicating root',
     *       'finding required')
99985 FORMAT ('    conflicts with the positive value NEQG,',I16,', pas',
     *       'sed to D02QWF.')
99984 FORMAT (' ** On entry, TOUT = T: T is ',1P,D13.5,'.')
99983 FORMAT (' ** On entry, TOUT = ',1P,D13.5,', but CRIT was set .TR',
     *       'UE. in the call')
99982 FORMAT ('    to D02QWF and integration cannot be attempted beyon',
     *       'd TCRIT,',1P,D13.5,'.')
99981 FORMAT (' ** The value of the argument T has been changed from ',
     *       1P,D13.5,' to')
99980 FORMAT ('    ',1P,D13.5,'. This is not permitted on a continuati',
     *       'on call.')
99979 FORMAT (' ** The input value of TOUT,',1P,D13.5,', indicates a c',
     *       'hange in the ')
99978 FORMAT ('    integration direction. This is not permitted on a c',
     *       'ontinuation call.')
99977 FORMAT (' ** The maximum number of steps has been attempted.')
99976 FORMAT (' ** The error tolerances are too stringent.')
99975 FORMAT (' ** ATOL(I) was set to 0.0E0 and now Y(I) is 0.0E0, for',
     *       ' some I.')
99974 FORMAT (' ** The problem appears to be stiff.')
99973 FORMAT (' ** Apparent convergence to a singular point of an even',
     *       't equation.')
99972 FORMAT (' ** Two successive errors detected at the current value',
     *       ' of T,',1P,D13.5,'.')
99971 FORMAT (' ** The input value of IREVCM,',I16,', is illegal.')
      END
