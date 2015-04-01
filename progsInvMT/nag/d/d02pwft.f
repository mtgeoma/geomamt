      SUBROUTINE D02PWF(TENDNU,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE RESET $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code D02PWF and how it is used in
C  conjunction with D02PDF to solve initial value problems, you should
C  study the document file RKSUITE.DOC carefully before attempting to
C  use the code. The following "Brief Reminder" is intended only to
C  remind you of the meaning, type, and size requirements of the
C  arguments.
C
C  The integration using D02PDF proceeds from TSTART in the direction of
C  TEND, and is now at TNOW.  To reset TEND to a new value TENDNU, just
C  call D02PWF with TENDNU as the argument.  You must continue
C  integrating in the same direction, so the sign of (TENDNU - TNOW)
C  must be the same as the sign of (TEND - TSTART). To change direction
C  you must restart by a call to D02PVF.
C
C  INPUT VARIABLE
C
C     TENDNU    - DOUBLE PRECISION
C                 The new value of TEND.
C                 Constraint: TENDNU and TNOW must be clearly
C                 distinguishable in the precision used.  The
C                 sign of (TENDNU - TNOW) must be the same as
C                 the sign of (TEND - TSTART).
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02PWF')
      LOGICAL           ASK
      INTEGER           MINUS1, MINUS2
      PARAMETER         (ASK=.TRUE.,MINUS1=-1,MINUS2=-2)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TENDNU
      INTEGER           IFAIL
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, CUBRMC, DIR, DWARF, EXPON, H, HOLD, HSTRT,
     *                  MCHEPS, RNDOFF, RS, RS1, RS2, RS3, RS4, SAFETY,
     *                  SQRRMC, STBRAD, T, TANANG, TINY, TND, TOLD,
     *                  TOLR, TOOSML, TSTRT
      INTEGER           FLSTP, LSTSTG, MAXTRY, NEQN, NFCN, NSEC, OKSTP,
     *                  ORDER, OUTCH, SVNFCN
      LOGICAL           FIRST, FSAL, LAST, UTASK
C     .. Arrays in Common ..
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  HMIN, TDIFF
      INTEGER           IER, NREC, STATE, STATE1
C     .. External Subroutines ..
      EXTERNAL          D02PDM, D02PDP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /ED02PD/TOOSML, COST, SAFETY, EXPON, STBRAD,
     *                  TANANG, RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG,
     *                  MAXTRY, NSEC, FSAL
      COMMON            /GD02PD/MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC,
     *                  TINY, OUTCH
      COMMON            /HD02PD/UTASK
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /AD02PD/, /BD02PD/, /ED02PD/, /GD02PD/,
     *                  /HD02PD/, /JD02PD/
C     .. Executable Statements ..
      IER = 0
      NREC = 0
C
C  Is it permissible to call D02PWF?
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
     *  ' ** You have called D02PWF after you specified to D02PVF that '
     *        ,
     *        ' ** you were going to use D02PCF. This is not permitted.'
            GO TO 20
         END IF
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *     ' ** You have not called D02PDF, so you cannot use D02PWF.'
         GO TO 20
      END IF
      IF (STATE.EQ.5 .OR. STATE.EQ.6) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(A,I1,A/A)')
     *     ' ** D02PDF has returned with IFAIL =  ', STATE, '.',
     *     ' ** You cannot call D02PWF in this circumstance.'
         GO TO 20
      END IF
C
C  Check value of TENDNU
C
      IF (DIR.GT.ZERO .AND. TENDNU.LE.T) THEN
         IER = 1
         NREC = 4
         WRITE (REC,FMT='(A/A,D13.5/A,D13.5,A/A)')
     *   ' ** Integration is proceeding in the positive direction. The '
     *     , ' ** current value for the independent variable is ', T,
     *     ' ** and you have set TENDNU = ', TENDNU,
     *     '.  TENDNU must be ', ' ** greater than T.'
      ELSE IF (DIR.LT.ZERO .AND. TENDNU.GE.T) THEN
         IER = 1
         NREC = 4
         WRITE (REC,FMT='(A/A,D13.5/A,D13.5,A/A)')
     *   ' ** Integration is proceeding in the negative direction. The '
     *     , ' ** current value for the independent variable is ', T,
     *     ' ** and you have set TENDNU = ', TENDNU,
     *     '.  TENDNU must be ', ' ** less than T.'
      ELSE
         HMIN = MAX(TINY,TOOSML*MAX(ABS(T),ABS(TENDNU)))
         TDIFF = ABS(TENDNU-T)
         IF (TDIFF.LT.HMIN) THEN
            IER = 1
            NREC = 4
            WRITE (REC,FMT='(A,D13.5,A/A,D13.5,A/A/A,D13.5,A)')
     *        ' ** The current value of the independent variable T is ',
     *        T, '.', ' ** The TENDNU you supplied has ABS(TENDNU-T) = '
     *        , TDIFF, '.',
     *     ' ** For the METHOD and the precision of the computer being '
     *        , ' ** used, this difference must be at least ', HMIN, '.'
         END IF
      END IF
      IF (IER.EQ.1) GO TO 20
C
      TND = TENDNU
      LAST = .FALSE.
C
   20 CONTINUE
C
      CALL D02PDP(IER,SRNAME,NREC,IFAIL)
C
      RETURN
      END
