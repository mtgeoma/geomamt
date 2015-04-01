      SUBROUTINE D02PYF(TOTFCN,STPCST,WASTE,STPSOK,HNEXT,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE STAT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code D02PYF and how it is used in
C  conjunction with the integrators D02PDF and D02PCF, you should study
C  the document file RKSUITE.DOC carefully before attempting to use the
C  code. The following "Brief Reminder" is intended only to remind you
C  of the meaning, type, and size requirements of the arguments.
C
C  D02PYF is called to obtain some details about the integration.
C
C  OUTPUT VARIABLES
C
C     TOTFCN    - INTEGER
C                 Total number of calls to F in the integration so far
C                 -- a measure of the cost of the integration.
C     STPCST    - INTEGER
C                 Cost of a typical step with this METHOD measured in
C                 calls to F.
C     WASTE     - DOUBLE PRECISION
C                 The number of attempted steps that failed to meet the
C                 local error requirement divided by the total number of
C                 steps attempted so far in the integration.
C     STPSOK    - INTEGER
C                 The number of accepted steps.
C     HNEXT     - DOUBLE PRECISION
C                 The step size the integrator plans to use for the next
C                 step.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02PYF')
      LOGICAL           ASK
      INTEGER           MINUS1, MINUS2
      PARAMETER         (ASK=.TRUE.,MINUS1=-1,MINUS2=-2)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HNEXT, WASTE
      INTEGER           IFAIL, STPCST, STPSOK, TOTFCN
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, EXPON, H, HOLD, LOCMAX, MAXERR, RS, RS1,
     *                  RS2, RS3, RS4, SAFETY, STBRAD, T, TANANG, TOLD,
     *                  TOOSML
      INTEGER           FLSTP, GNFCN, LSTSTG, MAXTRY, NFCN, NSEC, OKSTP,
     *                  ORDER, PRZERR, PRZERS, PRZSTG, PRZY, PRZYNU,
     *                  PRZYP, SVNFCN
      LOGICAL           ERASFL, ERASON, FIRST, FSAL, LAST, UTASK
C     .. Arrays in Common ..
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      INTEGER           IER, NREC, STATE
C     .. External Subroutines ..
      EXTERNAL          D02PDM, D02PDP
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Common blocks ..
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /ED02PD/TOOSML, COST, SAFETY, EXPON, STBRAD,
     *                  TANANG, RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG,
     *                  MAXTRY, NSEC, FSAL
      COMMON            /FD02PD/MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY,
     *                  PRZYP, PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
      COMMON            /HD02PD/UTASK
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /BD02PD/, /ED02PD/, /FD02PD/, /HD02PD/, /JD02PD/
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
C  Is it permissible to call D02PYF?
C
      CALL D02PDM(ASK,SRNAME,STATE)
      IF (STATE.EQ.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *   ' ** A catastrophic error has already been detected elsewhere.'
         GO TO 20
      END IF
      IF (STATE.EQ.MINUS2) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A/A)')
     *     ' ** You have already made one call to D02PYF after the',
     *    ' ** integrator returned with IFAIL set to 5 or 6. You cannot'
     *     , ' ** call D02PYF again.'
         GO TO 20
      END IF
      CALL D02PDM(ASK,'D02PDF',STATE)
      IF (STATE.EQ.MINUS1) THEN
         IER = 1
         NREC = 1
         IF (UTASK) THEN
            WRITE (REC,FMT='(A)')
     *       ' ** You have not called D02PCF, so you cannot use D02PYF.'
         ELSE
            WRITE (REC,FMT='(A)')
     *       ' ** You have not called D02PDF, so you cannot use D02PYF.'
         END IF
         GO TO 20
      END IF
C
C  Set flag so that the routine can only be called once after a hard
C  failure from the integrator.
      IF (STATE.EQ.5 .OR. STATE.EQ.6) IER = MINUS2
C
      TOTFCN = SVNFCN + NFCN
      IF (ERASON) TOTFCN = TOTFCN + GNFCN
      STPCST = COST
      STPSOK = OKSTP
      IF (OKSTP.LE.1) THEN
         WASTE = ZERO
      ELSE
         WASTE = DBLE(FLSTP)/DBLE(FLSTP+OKSTP)
      END IF
      HNEXT = H
C
   20 CONTINUE
C
      CALL D02PDP(IER,SRNAME,NREC,IFAIL)
C
      RETURN
      END
