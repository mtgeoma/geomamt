      SUBROUTINE D02QXF(NEQF,YP,TCURR,HLAST,HNEXT,ODLAST,ODNEXT,NSUCC,
     *                  NFAIL,TOLFAC,BADCMP,RWORK,LRWORK,IWORK,LIWORK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      INTEGER           SVNEQF, SVLRWK, SVLIWK, SVNSUC, SVNFAI, SVBADC,
     *                  SVKOLD, SVKORD, HOLDSV, HNXTSV, TLFCSV, TCURSV
      PARAMETER         (SVNEQF=1,SVLRWK=6,SVLIWK=7,SVNSUC=13,SVNFAI=14,
     *                  SVBADC=16,SVKOLD=18,SVKORD=19,HOLDSV=10,
     *                  HNXTSV=11,TLFCSV=12,TCURSV=13)
      INTEGER           ASK, NOTOK
      PARAMETER         (ASK=0,NOTOK=0)
      INTEGER           MAXREC
      PARAMETER         (MAXREC=6)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02QXF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HLAST, HNEXT, TCURR, TOLFAC
      INTEGER           BADCMP, IFAIL, LIWORK, LRWORK, NEQF, NFAIL,
     *                  NSUCC, ODLAST, ODNEXT
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), YP(NEQF)
      INTEGER           IWORK(LIWORK)
C     .. Local Scalars ..
      INTEGER           IER, K, NREC, STATE
C     .. Local Arrays ..
      CHARACTER*80      REC(MAXREC)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02QFY
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
         END IF
      END IF
C
      IF (IER.NE.1) THEN
         ODLAST = IWORK(SVKOLD)
         ODNEXT = IWORK(SVKORD)
         HLAST = RWORK(HOLDSV)
         HNEXT = RWORK(HNXTSV)
         TCURR = RWORK(TCURSV)
         TOLFAC = RWORK(TLFCSV)
         BADCMP = IWORK(SVBADC)
         NSUCC = IWORK(SVNSUC)
         NFAIL = IWORK(SVNFAI)
         DO 20 K = 1, NEQF
            YP(K) = RWORK(20+K)
   20    CONTINUE
      END IF
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** Neither of the integrator routines D02QFF/D02QGF  h',
     *       'as been called.')
99998 FORMAT (' ** The dimension of RWORK supplied to D02QXF is not th',
     *       'at supplied to D02QWF:')
99997 FORMAT ('    LRWORK in D02QWF was ',I16,', LRWORK in D02QXF is ',
     *       I16,'.')
99996 FORMAT (' ** The dimension of IWORK supplied to D02QXF is not th',
     *       'at supplied to D02QWF:')
99995 FORMAT ('    LIWORK in D02QWF was ',I16,', LIWORK in D02QXF is ',
     *       I16,'.')
99994 FORMAT (' ** The value of NEQF supplied to D02QXF is not that su',
     *       'pplied to D02QWF:')
99993 FORMAT ('    NEQF in D02QWF was ',I16,', NEQF in D02QXF is ',I16,
     *       '.')
      END
