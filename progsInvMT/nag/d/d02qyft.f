      SUBROUTINE D02QYF(NEQG,INDEX,TYPE,EVENTS,RESIDS,RWORK,LRWORK,
     *                  IWORK,LIWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      INTEGER           SVNEQF, SVNEQG, SVVCTL, SVLRWK, SVLIWK, SVROOT
      PARAMETER         (SVNEQF=1,SVNEQG=2,SVVCTL=9,SVLRWK=6,SVLIWK=7,
     *                  SVROOT=15)
      INTEGER           ASK, NOTOK
      PARAMETER         (ASK=0,NOTOK=0)
      INTEGER           MAXREC
      PARAMETER         (MAXREC=7)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02QYF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, INDEX, LIWORK, LRWORK, NEQG, TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  RESIDS(NEQG), RWORK(LRWORK)
      INTEGER           EVENTS(NEQG), IWORK(LIWORK)
C     .. Local Scalars ..
      INTEGER           IBASE, IER, K, NREC, STATE
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
         IF (IWORK(SVROOT).EQ.0) THEN
            IER = 1
            NREC = 1
            WRITE (REC(1),FMT=99998)
         END IF
         IF (IWORK(SVLRWK).NE.LRWORK) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99997)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99996) IWORK(SVLRWK), LRWORK
         END IF
         IF (IWORK(SVLIWK).NE.LIWORK) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99995)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99994) IWORK(SVLIWK), LIWORK
         END IF
         IF (IWORK(SVNEQG).NE.NEQG) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99993)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99992) IWORK(SVNEQG), NEQG
         END IF
      END IF
C
      IF (IER.NE.1) THEN
         INDEX = IWORK(11)
         TYPE = IWORK(12)
         DO 20 K = 1, NEQG
            EVENTS(K) = IWORK(20+K)
   20    CONTINUE
         IF (IWORK(SVVCTL).EQ.1) THEN
            IBASE = 20 + 21*IWORK(SVNEQF) + 2*IWORK(SVNEQF)
         ELSE
            IBASE = 20 + 21*IWORK(SVNEQF) + 2
         END IF
         DO 40 K = 1, NEQG
            RESIDS(K) = RWORK(IBASE+K)
   40    CONTINUE
      END IF
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** Neither of the integrator routines D02QFF/D02QGF ha',
     *       's been called.')
99998 FORMAT (' ** The integrator did not end at an event.')
99997 FORMAT (' ** The dimension of RWORK supplied to D02QYF is not th',
     *       'at supplied to D02QWF:')
99996 FORMAT ('    LRWORK in D02QWF was ',I16,', LRWORK in D02QYF is ',
     *       I16,'.')
99995 FORMAT (' ** The dimension of IWORK supplied to D02QYF is not th',
     *       'at supplied to D02QWF:')
99994 FORMAT ('    LIWORK in D02QWF was ',I16,', LIWORK in D02QYF is ',
     *       I16,'.')
99993 FORMAT (' ** The value of NEQG supplied to D02QYF is not that su',
     *       'pplied to D02QWF:')
99992 FORMAT ('    NEQG in D02QWF was ',I16,', NEQG in D02QYF is ',I16,
     *       '.')
      END
