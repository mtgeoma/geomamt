      SUBROUTINE D02CGF(T,TEND,NEQ,Y,TOL,HMAX,M,VAL,FCN,RWORK,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02CGF')
      DOUBLE PRECISION  RZERO, ONE, HALF, TWO
      PARAMETER         (RZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TWO=2.0D0)
      INTEGER           MXREC
      PARAMETER         (MXREC=4)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HMAX, T, TEND, TOL, VAL
      INTEGER           IFAIL, M, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(28+21*NEQ), Y(NEQ)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
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
C     .. Arrays in Common ..
      DOUBLE PRECISION  ALPHA(12), BETA(12), GI(11), GRWB(13), PSI(12),
     *                  SIG(13), V(12), W(12)
      INTEGER           IV(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIR, GWANT, HINC, HMAX1, TLAST, TOUT, TSTOP,
     *                  TWANT
      INTEGER           BADCMP, I, IDID, IERR, IGNEW, IGOLD, IGP, IGSC,
     *                  IGV, IINDXG, IMMR, INEEDG, IP, IPHI, IPHI3N,
     *                  IPHI4N, IPHI5N, IPRO, IREVCM, IROD, ISEE, ITGV,
     *                  ITK, ITLB, ITRB, IWT, IYP, IYPOUT, IYY, KWANT,
     *                  MAXNUM, NEQG, NFAIL, NREC, NSUCC
      LOGICAL           TSTINC
C     .. Local Arrays ..
      DOUBLE PRECISION  ATOL(1), RTOL(1)
      INTEGER           IWORK(21)
      CHARACTER*80      REC(MXREC)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02QFX
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
     *                  GSTOP, PGSTOP, ROOT, ROOTS, NEEDG, DISCOP,
     *                  NEWGEQ, SEARCH, PSERCH
C     .. Save statement ..
      SAVE              /AD02QF/, /BD02QF/
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
      IF (M.LT.1) THEN
         IERR = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99996) M
      ELSE IF (M.GT.NEQ) THEN
         IERR = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99995) M, NEQ
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
      RTOL(1) = TOL
      ATOL(1) = TOL
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
      GSTOP = .TRUE.
      SEARCH = .FALSE.
      NEWGEQ = .TRUE.
      NEQG = 1
      ISEE = 0
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
      TSTAR = T
      TLAST = T
      HMAX1 = RZERO
      MAXNUM = 1000
      NSUCC = 0
      NFAIL = 0
      HINC = ABS(HMAX)
      DIR = SIGN(ONE,TEND-T)
      TSTINC = HINC .NE. RZERO .AND. HINC .LE. ABS(TEND-T)
      IF (TSTINC) THEN
         HINC = SIGN(HINC,DIR)
         IF (DIR*(T+TWO*HINC-TEND).GE.RZERO) THEN
            TOUT = T + HALF*(TEND-T)
         ELSE
            TOUT = T + HINC
         END IF
      ELSE
         TOUT = TEND
      END IF
C
      IREVCM = 0
   20 CONTINUE
      CALL D02QFX(IREVCM,TWANT,KWANT,GWANT,NEQ,T,Y,TOUT,RTOL,ATOL,IDID,
     *            RWORK(IYPOUT),RWORK(IYP),RWORK(IYY),RWORK(IWT),
     *            RWORK(IP),RWORK(IPHI),RWORK(IGOLD),RWORK(IGNEW),
     *            RWORK(ITGV),RWORK(IGV),RWORK(IGP),RWORK(ITK),
     *            RWORK(ITLB),RWORK(ITRB),RWORK(IPRO),RWORK(IROD),
     *            RWORK(1),RWORK(11),RWORK(12),RWORK(13),HMAX1,MAXNUM,
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
         RWORK(IGOLD+I-1) = Y(M) - VAL
  160 CONTINUE
      GO TO 20
  180 CONTINUE
      GWANT = Y(M) - VAL
      GO TO 20
C
  200 CONTINUE
      IF (IDID.EQ.2 .OR. IDID.EQ.3) THEN
         IF (TSTINC .AND. TOUT.NE.TEND) THEN
            IF (DIR*(TOUT+HINC-TEND).GE.RZERO) THEN
               TOUT = TEND
            ELSE IF (DIR*(TOUT+TWO*HINC-TEND).GE.RZERO) THEN
               TOUT = TOUT + HALF*(TEND-TOUT)
            ELSE
               TOUT = TOUT + HINC
            END IF
            GO TO 20
         ELSE
            IERR = 4
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99994)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99993)
         END IF
      ELSE IF (IDID.EQ.-1 .OR. IDID.EQ.-4) THEN
         IBEGIN = 1
         GO TO 20
      ELSE IF (IDID.EQ.-2) THEN
         IF (T.EQ.TLAST) THEN
            IERR = 3
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99992)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99991) TOL
         ELSE
            IERR = 2
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99990) T
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99989) TOL
         END IF
      ELSE IF (IDID.LT.4 .OR. IDID.GT.6) THEN
         IERR = 5
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99988) IDID
      END IF
  220 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** N .lt. 1. N = ',I16,'.')
99998 FORMAT (' ** TOL .le. 0.0. TOL = ',1P,D13.5,'.')
99997 FORMAT (' ** X .eq. XEND. X = ',1P,D13.5,'.')
99996 FORMAT (' ** M .lt. 1. M = ',I16,'.')
99995 FORMAT (' ** M .gt. N. M = ',I16,', N = ',I16,'.')
99994 FORMAT (' ** No change in sign of the function Y(M) - VAL was de',
     *       'tected ')
99993 FORMAT ('    in the integration range.')
99992 FORMAT (' ** No integration steps have been taken. Progress not ',
     *       'possible with the ')
99991 FORMAT ('    input value of TOL,',1P,D13.5,'.')
99990 FORMAT (' ** Integration successful as far as T = ',1P,D13.5,', ',
     *       'but further progress')
99989 FORMAT ('    not possible with the input value of TOL,',1P,D13.5,
     *       '.')
99988 FORMAT (' ** Impossible error - internal variable IDID = ',I16,
     *       '.')
      END
