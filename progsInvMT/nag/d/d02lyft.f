      SUBROUTINE D02LYF(NEQ,HNEXT,HUSED,HSTART,NSUCC,NFAIL,NATT,THRES,
     *                  THRESP,RWORK,LRWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     THIS IS THE DIAGNOSTIC ROUTINE ASSOCIATED WITH THE NYSTROM
C     SOLVER D02LAF.  THE VALUES OF THE PARAMETERS IN THE SUBROUTINE
C     CALL ARE TAKEN OUT OF RWORK.
C
C     SEE COMMENT IN D02LXF ABOUT STORAGE REDUCTION FOR THRES AND
C     THRESP.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02LYF')
      DOUBLE PRECISION  XINC
      INTEGER           NOVHD, ASK, NOTOK
      PARAMETER         (XINC=0.2D0,NOVHD=16,ASK=0,NOTOK=0)
      INTEGER           SVINIH, SVHUSD, SVHNXT, SVOKST, SVFLST, SVATST,
     *                  SVNEQ, SVLRWK
      PARAMETER         (SVINIH=1,SVHUSD=2,SVHNXT=3,SVOKST=4,SVFLST=5,
     *                  SVATST=11,SVNEQ=15,SVLRWK=16)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HNEXT, HSTART, HUSED
      INTEGER           IFAIL, LRWORK, NATT, NEQ, NFAIL, NSUCC
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), THRES(NEQ), THRESP(NEQ)
C     .. Local Scalars ..
      INTEGER           IER, K, NREC, NTHR, NTHRP, STATE
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02LAX
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
      CALL D02LAX(STATE,ASK)
C
      IF (STATE.EQ.NOTOK) THEN
         IER = 1
         NREC = 1
         WRITE (REC(1),FMT=99999)
      ELSE
         IF (NEQ.NE.INT(RWORK(SVNEQ))) THEN
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99998)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99997) INT(RWORK(SVNEQ)), NEQ
            IER = 1
         END IF
         IF (LRWORK.NE.INT(RWORK(SVLRWK))) THEN
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99996)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99995) INT(RWORK(SVLRWK)), LRWORK
            IER = 1
         END IF
      END IF
C
      NTHR = NOVHD + 4*NEQ
      NTHRP = NOVHD + 5*NEQ
C
      IF (IER.NE.1) THEN
         HSTART = RWORK(SVINIH)
         HUSED = RWORK(SVHUSD)
         HNEXT = RWORK(SVHNXT)
         NSUCC = INT(RWORK(SVOKST)+XINC)
         NFAIL = INT(RWORK(SVFLST)+XINC)
         NATT = INT(RWORK(SVATST)+XINC)
         DO 20 K = 1, NEQ
            THRES(K) = RWORK(NTHR+K)
            THRESP(K) = RWORK(NTHRP+K)
   20    CONTINUE
      END IF
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** The integrator routine D02LAF has not been called.')
99998 FORMAT (' ** The value of NEQ supplied to D02LYF is not that sup',
     *       'plied to D02LXF:')
99997 FORMAT ('    NEQ in D02LXF was ',I16,', NEQ in D02LYF is ',I16,
     *       '.')
99996 FORMAT (' ** The dimension of RWORK supplied to D02LYF is not th',
     *       'at supplied to D02LXF:')
99995 FORMAT ('    LRWORK in D02LXF was ',I16,', LRWORK in D02LYF is ',
     *       I16,'.')
      END
