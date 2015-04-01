      SUBROUTINE D02PZF(RMSERR,ERRMAX,TERRMX,WORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE GLBERR $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code D02PZF and how it is used in
C  conjunction with D02PCF and D02PDF to solve initial value problems,
C  you should study the document file RKSUITE.DOC carefully before
C  attempting to use the code.  The following "Brief Reminder" is
C  intended only to remind you of the meaning, type, and size
C  requirements of the arguments.
C
C  If ERRASS was set .TRUE. in the call to D02PVF, then after any call
C  to D02PCF or D02PDF to advance the integration to TNOW or TWANT, the
C  subroutine D02PZF may be called to obtain an assessment of the true
C  error of the integration. At each step and for each solution
C  component Y(L), a more accurate "true" solution YT(L), an average
C  magnitude "size(L)" of its size, and its error
C  abs(Y(L) - YT(L))/max("size(L)",THRES(L)) are computed. The
C  assessment returned in RMSERR(L) is the RMS (root-mean-square)
C  average of the error in the Lth solution component over all steps of
C  the integration from TSTART through TNOW.
C
C  OUTPUT VARIABLES
C
C     RMSERR(*) - DOUBLE PRECISION array of length NEQ
C                 RMSERR(L) approximates the RMS average of the true
C                 error of the numerical solution for the Ith solution
C                 component, L = 1,2,...,NEQ.  The average is taken over
C                 all steps from TSTART to TNOW.
C     ERRMAX    - DOUBLE PRECISION
C                 The maximum (approximate) true error taken over all
C                 solution components and all steps from TSTART to TNOW.
C     TERRMX    - DOUBLE PRECISION
C                 First value of the independent variable where the
C                 (approximate) true error attains the maximum value
C                 ERRMAX.
C
C  WORKSPACE
C
C     WORK(*)   - DOUBLE PRECISION array as used in D02PVF and D02PCF or
C                 D02PDF
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02PZF')
      LOGICAL           ASK
      PARAMETER         (ASK=.TRUE.)
      INTEGER           MINUS1, MINUS2
      PARAMETER         (MINUS1=-1,MINUS2=-2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERRMAX, TERRMX
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  RMSERR(*), WORK(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DIR, H, HOLD, HSTRT, LOCMAX, MAXERR, T, TND,
     *                  TOLD, TOLR, TSTRT
      INTEGER           FLSTP, GNFCN, NEQN, NFCN, OKSTP, PRZERR, PRZERS,
     *                  PRZSTG, PRZY, PRZYNU, PRZYP, SVNFCN
      LOGICAL           ERASFL, ERASON, FIRST, LAST, UTASK
C     .. Arrays in Common ..
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      INTEGER           IER, L, NREC, STATE
C     .. External Subroutines ..
      EXTERNAL          D02PDM, D02PDP
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /FD02PD/MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY,
     *                  PRZYP, PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
      COMMON            /HD02PD/UTASK
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /AD02PD/, /BD02PD/, /FD02PD/, /HD02PD/, /JD02PD/
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
C  Is it permissible to call D02PZF?
C
      CALL D02PDM(ASK,SRNAME,STATE)
      IF (STATE.EQ.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *   ' ** A catastrophic error has already been detected elsewhere.'
         GO TO 40
      END IF
      IF (STATE.EQ.MINUS2) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A/A)')
     *     ' ** You have already made one call to D02PZF after the',
     *    ' ** integrator returned with IFAIL set to 5 or 6. You cannot'
     *     , ' ** call D02PZF again.'
         GO TO 40
      END IF
      CALL D02PDM(ASK,'D02PDF',STATE)
      IF (STATE.EQ.MINUS1) THEN
         IER = 1
         NREC = 1
         IF (UTASK) THEN
            WRITE (REC,FMT='(A)')
     *  ' ** You have not yet called D02PCF, so you cannot call D02PZF.'
         ELSE
            WRITE (REC,FMT='(A)')
     *  ' ** You have not yet called D02PDF, so you cannot call D02PZF.'
         END IF
         GO TO 40
      END IF
C
C  Set flag so that the routine can only be called once after a hard
C  failure from the integrator.
      IF (STATE.EQ.5 .OR. STATE.EQ.6) IER = MINUS2
C
C  Check that ERRASS was set properly for error assessment in D02PVF.
C
      IF ( .NOT. ERASON) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A/A)')
     *     ' ** No error assessment is available since you did not ',
     *     ' ** ask for it in your call to the routine D02PVF.',
     *     ' ** Check your program carefully.'
         GO TO 40
      END IF
C
C  Check to see if the integrator has not actually taken a step.
C
      IF (OKSTP.EQ.0) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A/A)')
     *     ' ** The integrator has not actually taken any successful ',
     *     ' ** steps.  You cannot call D02PZF in this circumstance. ',
     *     ' ** Check your program carefully.'
         GO TO 40
      END IF
C
C  Compute RMS error and set output variables.
C
      ERRMAX = MAXERR
      TERRMX = LOCMAX
      DO 20 L = 1, NEQN
         RMSERR(L) = SQRT(WORK(PRZERR-1+L)/OKSTP)
   20 CONTINUE
C
   40 CONTINUE
C
      CALL D02PDP(IER,SRNAME,NREC,IFAIL)
C
      RETURN
      END
