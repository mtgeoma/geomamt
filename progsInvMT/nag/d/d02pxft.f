      SUBROUTINE D02PXF(TWANT,REQEST,NWANT,YWANT,YPWANT,F,WORK,WRKINT,
     *                  LENINT,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE INTRP $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code D02PXF and how it is used in
C  conjunction with D02PDF to solve initial value problems, you should
C  study the document file RKSUITE.DOC carefully before attempting to
C  use the code. The following "Brief Reminder" is intended only to
C  remind you of the meaning, type, and size requirements of the
C  arguments.
C
C  When integrating with METHOD = 1 or 2, answers may be obtained
C  inexpensively between steps by interpolation. D02PXF is called after
C  a step by D02PDF from a previous value of T, called TOLD below, to
C  the current value of T to get an answer at TWANT. You can specify any
C  value of TWANT you wish, but specifying a value outside the interval
C  [TOLD,T] is likely to yield answers with unsatisfactory accuracy.
C
C  INPUT VARIABLE
C
C     TWANT     - DOUBLE PRECISION
C                 The value of the independent variable where a solution
C                 is desired.
C
C  The interpolant is to be evaluated at TWANT to approximate the
C  solution and/or its first derivative there.  There are three
C  cases:
C
C  INPUT VARIABLE
C
C     REQEST    - CHARACTER*(*)
C                 Only the first character of REQEST is significant.
C                 REQEST(1:1) = `S' or `s'- compute approximate
C                             = `D' or `d'- compute approximate first
C                                           `D'erivative of the solution
C                                           only.
C                             = `B' or `b'- compute `B'oth approximate
C                                           solution and first
C                                           derivative.
C                 Constraint: REQEST(1:1) must be `S',`s',`D',`d',`B' or
C                             `b'.
C
C  If you intend to interpolate at many points, you should arrange for
C  the the interesting components to be the first ones because the code
C  approximates only the first NWANT components.
C
C  INPUT VARIABLE
C
C     NWANT     - INTEGER
C                 Only the first NWANT components of the answer are to
C                 be computed.
C                 Constraint:  NEQ >= NWANT >= 1
C
C  OUTPUT VARIABLES
C
C     YWANT(*)  - DOUBLE PRECISION array of length NWANT
C                 Approximation to the first NWANT components of the
C                 true solution at TWANT when REQESTed.
C     YPWANT(*) - DOUBLE PRECISION array of length NWANT
C                 Approximation to the first NWANT components of the
C                 first derivative of the true solution at TWANT when
C                 REQESTed.
C
C  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE PROGRAM CALLING D02PXF:
C
C     F         - name of the subroutine for evaluating the differential
C                 equations as provided to D02PDF.
C
C  WORKSPACE
C
C     WORK(*)   - DOUBLE PRECISION array as used in D02PVF and D02PDF
C                 Do not alter the contents of this array.
C
C     WRKINT(*) - DOUBLE PRECISION array of length LENINT
C                 Do not alter the contents of this array.
C
C     LENINT    - INTEGER
C                 Length of WRKINT. If
C                 METHOD = 1, LENINT must be at least 1
C                        = 2, LENINT must be at least
C                          NEQ+MAX(NEQ,5*NWANT)
C                        = 3--CANNOT BE USED WITH THIS SUBROUTINE
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02PXF')
      LOGICAL           ASK
      INTEGER           PLUS1, MINUS1, MINUS2
      PARAMETER         (ASK=.TRUE.,PLUS1=1,MINUS1=-1,MINUS2=-2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TWANT
      INTEGER           IFAIL, LENINT, NWANT
      CHARACTER*(*)     REQEST
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*), WRKINT(LENINT), YPWANT(*), YWANT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  DIR, HSTRT, TND, TOLR, TSTRT
      INTEGER           LNINTP, METHD, MINTP, NEQN, NSTAGE, PRERST,
     *                  PRINTP, PRSCR, PRSTGS, PRTHRS, PRWT, PRY,
     *                  PRYOLD, PRYP
      LOGICAL           INTP, UTASK
C     .. Arrays in Common ..
      DOUBLE PRECISION  A(13,13), B(13), BHAT(13), C(13), E(7), R(11,6)
      INTEGER           PTR(13)
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      INTEGER           ICHK, IER, NREC, NWNTSV, STARTP, STATE, STATE1
      LOGICAL           ININTP, LEGALR
      CHARACTER         REQST1
C     .. External Subroutines ..
      EXTERNAL          D02PDM, D02PDP, D02PXY, D02PXZ
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /CD02PD/PRTHRS, PRERST, PRWT, PRYOLD, PRSCR,
     *                  PRY, PRYP, PRSTGS, PRINTP, LNINTP
      COMMON            /DD02PD/A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
     *                  MINTP, INTP
      COMMON            /HD02PD/UTASK
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /AD02PD/, /CD02PD/, /DD02PD/, /HD02PD/,
     *                  /JD02PD/, NWNTSV, ININTP, STARTP
C     .. Data statements ..
      DATA              NWNTSV/MINUS1/
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
C  Is it permissible to call D02PXF?
C
      CALL D02PDM(ASK,'D02PDF',STATE)
      IF (STATE.EQ.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *   ' ** A catastrophic error has already been detected elsewhere.'
         GO TO 20
      END IF
      IF (UTASK) THEN
         CALL D02PDM(ASK,'D02PCF',STATE1)
         IF (STATE1.NE.MINUS2) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(A/A)')
     *       ' ** You have called D02PXF after you specified to D02PVF '
     *        ,
     *   ' ** that you were going to use D02PCF. This is not permitted.'
            GO TO 20
         END IF
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *     ' ** You have not called D02PDF, so you cannot use D02PXF.'
         GO TO 20
      END IF
      IF (STATE.GE.PLUS1) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(A/A)')
     *     ' ** D02PDF has returned with an IFAIL set greater than 0.',
     *     ' ** You cannot call D02PXF in this circumstance.'
         GO TO 20
      END IF
C
C  Check input
C
      REQST1 = REQEST(1:1)
      LEGALR = REQST1 .EQ. 'S' .OR. REQST1 .EQ. 's' .OR. REQST1 .EQ.
     *         'D' .OR. REQST1 .EQ. 'd' .OR. REQST1 .EQ. 'B' .OR.
     *         REQST1 .EQ. 'b'
      IF ( .NOT. LEGALR) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A,A,A/A)')
     *     ' ** You have set the first character of ',
     *     ' ** REQEST to be ''', REQST1, '''. It must be one of ',
     *     ' ** ''S'',''s'',''D'',''d'',''B'' or ''b''.'
         GO TO 20
      END IF
C
      IF (NWANT.GT.NEQN) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A,I6,A/A,I6,A/A)')
     *     ' ** You have specified the value of NWANT to be ', NWANT,
     *     '. This', ' ** is greater than ', NEQN,
     *     ', which is the number of equations ',
     *     ' ** in the system being integrated.'
         GO TO 20
      ELSE IF (NWANT.LT.1) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A,I6,A/A/A)')
     *     ' ** You have specified the value of NWANT to be ', NWANT,
     *     ', but ',
     *     ' ** this is less than 1. You cannot interpolate a zero or ',
     *     ' ** negative number of components.'
         GO TO 20
      END IF
C
      IF (METHD.EQ.1) THEN
         IF (LENINT.LT.1) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(A,I6,A/A)')
     *        ' ** You have specified LENINT to be ', LENINT, '.',
     *        ' ** This is too small. LENINT must be at least 1.'
            GO TO 20
         END IF
         STARTP = 1
      ELSE IF (METHD.EQ.2) THEN
         ICHK = NEQN + 5*NWANT
         IF (LENINT.LT.ICHK) THEN
            IER = 1
            NREC = 3
            WRITE (REC,FMT='(A,I6,A/A/A,I6,A)')
     *        ' ** You have specified LENINT to be ', LENINT,
     *        '. This is too',
     *        ' ** small. LENINT must be at least NEQ+5*NWANT which is',
     *        ' ** ', ICHK, '.'
            GO TO 20
         END IF
         STARTP = NEQN + 1
      ELSE IF (METHD.EQ.3) THEN
         IER = 1
         NREC = 5
         WRITE (REC,FMT='(A,A/A/A/A/A)')
     *    ' ** You have been using D02PDF with METHOD = 3 to integrate '
     *     , 'your',
     *  ' ** equations. You have just called D02PXF, but interpolation '
     *     ,
     *  ' ** is not available for this METHOD. Either use METHOD = 2,  '
     *     ,
     * ' ** for which interpolation is available, or use D02PWF to make'
     *     ,
     *    ' ** D02PDF step exactly to the points where you want output.'
         GO TO 20
      END IF
C
C  Has the interpolant been initialised for this step?
C
      CALL D02PDM(ASK,SRNAME,STATE)
      ININTP = STATE .NE. MINUS2
C
C  Some initialization must be done before interpolation is possible.
C  To reduce the overhead, the interpolating polynomial is formed for
C  the first NWANT components.  In the unusual circumstance that NWANT
C  is changed while still interpolating within the span of the current
C  step, the scheme must be reinitialized to accomodate the additional
C  components.
C
      IF ( .NOT. ININTP .OR. NWANT.NE.NWNTSV) THEN
C
C  At present the derivative of the solution at the previous step,
C  YPOLD(*), is stored in the scratch array area starting at PRSCR. In
C  the case of METHD = 1 we can overwrite the stages.
C
         IF (METHD.EQ.1) THEN
            CALL D02PXY(F,NEQN,NWANT,WORK(PRY),WORK(PRYP),WORK(PRYOLD),
     *                  WORK(PRSCR),WORK(PRSTGS), .NOT. ININTP,
     *                  WORK(PRSTGS),WORK(PRSTGS),WORK(PRSTGS))
         ELSE
            CALL D02PXY(F,NEQN,NWANT,WORK(PRY),WORK(PRYP),WORK(PRYOLD),
     *                  WORK(PRSCR),WORK(PRSTGS), .NOT. ININTP,WRKINT,
     *                  WORK(PRWT),WRKINT(STARTP))
         END IF
C
C  Set markers to show that interpolation has been initialized for
C  NWANT components.
         NWNTSV = NWANT
      END IF
C
C  The actual evaluation of the interpolating polynomial and/or its
C  first derivative is done in D02PXZ.
C
      IF (METHD.EQ.1) THEN
         CALL D02PXZ(WORK(PRY),WORK(PRYP),WORK(PRSTGS),TWANT,REQEST,
     *               NWANT,YWANT,YPWANT)
      ELSE
         CALL D02PXZ(WORK(PRY),WORK(PRYP),WRKINT(STARTP),TWANT,REQEST,
     *               NWANT,YWANT,YPWANT)
      END IF
C
   20 CONTINUE
C
      CALL D02PDP(IER,SRNAME,NREC,IFAIL)
C
      RETURN
      END
