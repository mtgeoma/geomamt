      SUBROUTINE D02QZF(NEQF,TWANT,NWANT,YWANT,YPWANT,RWORK,LRWORK,
     *                  IWORK,LIWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      INTEGER           SVNEQF, SVLRWK, SVLIWK, SVNSUC, HOLDSV, TCURSV
      PARAMETER         (SVNEQF=1,SVLRWK=6,SVLIWK=7,SVNSUC=13,HOLDSV=10,
     *                  TCURSV=13)
      INTEGER           ASK, NOTOK
      PARAMETER         (ASK=0,NOTOK=0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      INTEGER           MAXREC
      PARAMETER         (MAXREC=6)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02QZF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TWANT
      INTEGER           IFAIL, LIWORK, LRWORK, NEQF, NWANT
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), YPWANT(NWANT), YWANT(NWANT)
      INTEGER           IWORK(LIWORK)
C     .. Local Scalars ..
      INTEGER           IER, IP, IPHI, IYY, NREC, STATE
      LOGICAL           INTERP
C     .. Local Arrays ..
      CHARACTER*80      REC(MAXREC)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02QFY, D02QZZ
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
      CALL D02QFY(STATE,ASK)
C
      IF (STATE.EQ.NOTOK) THEN
         IER = 1
         NREC = 1
         WRITE (REC(1),FMT=99999)
      ELSE
         IF (IWORK(SVLRWK).NE.LRWORK) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99998)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99997) IWORK(SVLRWK), LRWORK
         END IF
         IF (IWORK(SVLIWK).NE.LIWORK) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99996)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99995) IWORK(SVLIWK), LIWORK
         END IF
         IF (IWORK(SVNEQF).NE.NEQF) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99994)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99993) IWORK(SVNEQF), NEQF
         ELSE IF (NWANT.GT.NEQF) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99990) NWANT, NEQF
         ELSE IF (NWANT.LT.1) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99989) NWANT
         END IF
      END IF
      IF (IER.EQ.1) GO TO 20
C
      IF (IWORK(SVNSUC).EQ.0) THEN
         IER = 1
         NREC = 1
         WRITE (REC(1),FMT=99992) TWANT
         GO TO 20
      END IF
C
      IF (SIGN(ONE,RWORK(HOLDSV)).EQ.ONE) THEN
         INTERP = RWORK(TCURSV) - RWORK(HOLDSV)
     *            .LE. TWANT .AND. TWANT .LE. RWORK(TCURSV)
      ELSE
         INTERP = RWORK(TCURSV) - RWORK(HOLDSV)
     *            .GE. TWANT .AND. TWANT .GE. RWORK(TCURSV)
      END IF
      IF ( .NOT. INTERP) THEN
         IER = 2
         NREC = 1
         WRITE (REC(1),FMT=99991) TWANT
      END IF
C
      IYY = 2*NEQF + 21
      IP = 2*NEQF + IYY
      IPHI = NEQF + IP
      CALL D02QZZ(TWANT,YWANT,YPWANT,NWANT,NEQF,RWORK(TCURSV),RWORK(IYY)
     *            ,RWORK(IP),RWORK(IPHI))
C
   20 CONTINUE
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** Neither of the integrator routines D02QFF/D02QGF ha',
     *       's been called.')
99998 FORMAT (' ** The dimension of RWORK supplied to D02QZF is not th',
     *       'at supplied to D02QWF:')
99997 FORMAT ('    LRWORK in D02QWF was ',I16,', LRWORK in D02QZF is ',
     *       I16,'.')
99996 FORMAT (' ** The dimension of IWORK supplied to D02QZF is not th',
     *       'at supplied to D02QWF:')
99995 FORMAT ('    LIWORK in D02QWF was ',I16,', LIWORK in D02QZF is ',
     *       I16,'.')
99994 FORMAT (' ** The value of NEQF supplied to D02QZF is not that su',
     *       'pplied to D02QWF:')
99993 FORMAT ('    NEQF in D02QWF was ',I16,', NEQF in D02QZF is ',I16,
     *       '.')
99992 FORMAT (' ** No successful steps, interpolation impossible at TW',
     *       'ANT = ',1P,D13.5,'.')
99991 FORMAT (' ** Extrapolation at TWANT = ',1P,D13.5,'.')
99990 FORMAT (' ** NWANT is greater than NEQF. NWANT =',I16,', NEQF =',
     *       I16,'.')
99989 FORMAT (' ** NWANT is less than 1. NWANT = ',I16,'.')
      END
