      SUBROUTINE D02QWF(STATEF,NEQF,VECTOL,ATOL,LATOL,RTOL,LRTOL,ONESTP,
     *                  CRIT,TCRIT,HMAX,MAXSTP,NEQG,ALTERG,SOPHST,RWORK,
     *                  LRWORK,IWORK,LIWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      CHARACTER         START, RSTART, CONTIN
      PARAMETER         (START='S',RSTART='R',CONTIN='C')
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02QWF')
      INTEGER           SVNEQF, SVNEQG, SVALTG, SV1STP, SVSPHS, SVLRWK,
     *                  SVLIWK, SVCRIT, SVVCTL, SVSTAT, SVMXST, CRITSV,
     *                  HMAXSV
      PARAMETER         (SVNEQF=1,SVNEQG=2,SVALTG=3,SV1STP=4,SVSPHS=5,
     *                  SVLRWK=6,SVLIWK=7,SVCRIT=8,SVVCTL=9,SVSTAT=10,
     *                  SVMXST=17,CRITSV=1,HMAXSV=2)
      INTEGER           MAXREC
      PARAMETER         (MAXREC=6)
      INTEGER           SET
      PARAMETER         (SET=1)
      DOUBLE PRECISION  ZERO, FOUR
      PARAMETER         (ZERO=0.0D0,FOUR=4.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HMAX, TCRIT
      INTEGER           IFAIL, LATOL, LIWORK, LRTOL, LRWORK, MAXSTP,
     *                  NEQF, NEQG
      LOGICAL           ALTERG, CRIT, ONESTP, SOPHST, VECTOL
      CHARACTER*1       STATEF
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(LATOL), RTOL(LRTOL), RWORK(LRWORK)
      INTEGER           IWORK(LIWORK)
C     .. Local Scalars ..
      DOUBLE PRECISION  FOURU
      INTEGER           ERSTAT, I, IBATL, IBRTL, IER, JI1, JR1, JR2,
     *                  LIWREQ, LRWREQ, NREC
      LOGICAL           FIRST
      CHARACTER*1       STATE
C     .. Local Arrays ..
      CHARACTER*80      REC(MAXREC)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02QWZ, E04UDU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Save statement ..
      SAVE              FIRST
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
      STATE = STATEF
      CALL E04UDU(STATE)
C
C     ON THE VERY FIRST CALL CHECK THAT STATEF = 'S'
C
      IF (FIRST) THEN
         IF (STATE.NE.START) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99999) STATE
            GO TO 100
         END IF
         FIRST = .FALSE.
      END IF
C
C     DATA CHECKS FOR A START
C
      IF (STATE.EQ.START) THEN
         IWORK(SVSTAT) = 0
C
C        CHECK FOR VALID NEQF
C
         IF (NEQF.LT.1) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99998) NEQF
            GO TO 100
         ELSE
            IWORK(SVNEQF) = NEQF
         END IF
C
C        CHECK IF ROOTFINDING CAPABILITY SPECIFIED
C
         IWORK(SVALTG) = 0
         IF (NEQG.LE.0) THEN
            IWORK(SVNEQG) = 0
         ELSE
            IWORK(SVNEQG) = NEQG
C
C           COPY VALUE OF SOPHST
C
            IF (SOPHST) THEN
               IWORK(SVSPHS) = 1
            ELSE
               IWORK(SVSPHS) = 0
            END IF
         END IF
C
C        COPY VALUE OF VECTOL AND CHECK FOR ATOL,RTOL SUFFICIENT LENGTH
C
         IF (VECTOL) THEN
            IWORK(SVVCTL) = 1
            IF (LRTOL.LT.NEQF) THEN
               IER = 1
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99997)
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99993) LRTOL, NEQF
            END IF
            IF (LATOL.LT.NEQF) THEN
               IER = 1
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99996)
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99992) LRTOL, NEQF
            END IF
         ELSE
            IWORK(SVVCTL) = 0
            IF (LRTOL.LT.1) THEN
               IER = 1
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99995)
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99991) LRTOL
            END IF
            IF (LATOL.LT.1) THEN
               IER = 1
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99994)
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99990) LATOL
            END IF
         END IF
C
C        DATA CHECKS FOR A RESTART OR CONTINUATION
C
      ELSE IF (STATE.EQ.RSTART .OR. STATE.EQ.CONTIN) THEN
         IF (STATE.EQ.RSTART) THEN
            IWORK(SVSTAT) = -1
         ELSE
            IWORK(SVSTAT) = 1
         END IF
C
C        CHECK NEQF UNALTERED
C
         IF (NEQF.NE.IWORK(SVNEQF)) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99989) STATE
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99988) NEQF, IWORK(SVNEQF)
         END IF
C
C        CHECK VECTOL UNALTERED
C
         IF (VECTOL .AND. IWORK(SVVCTL).EQ.0) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99987) STATE
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99986)
         ELSE IF ( .NOT. VECTOL .AND. IWORK(SVVCTL).EQ.1) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99987) STATE
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99985)
         END IF
C
C        CHECK IF ROOTFINDING CAPABILITY SPECIFIED
C
         IF (NEQG.LE.0) THEN
            IWORK(SVNEQG) = 0
         ELSE
C
C           COPY VALUE OF SOPHST
C
            IF (SOPHST) THEN
               IWORK(SVSPHS) = 1
            ELSE
               IWORK(SVSPHS) = 0
            END IF
C
C           CHECK FOR REDEFINITION OF EVENT FUNCTIONS
C
            IF (ALTERG) THEN
               IWORK(SVALTG) = 1
               IWORK(SVNEQG) = NEQG
            ELSE
               IWORK(SVALTG) = 0
               IF (NEQG.NE.IWORK(SVNEQG)) THEN
                  IER = 1
                  NREC = NREC + 1
                  WRITE (REC(NREC),FMT=99984) STATE
                  NREC = NREC + 1
                  WRITE (REC(NREC),FMT=99983) NEQG, IWORK(SVNEQG)
               END IF
            END IF
         END IF
      ELSE
C
C        STATEF INCORRECTLY SET
C
         IER = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99974) STATE
      END IF
      IF (IER.EQ.1) GO TO 100
C
C     CHECK WORKSPACE SUFFICIENT
C
      IF (NEQG.LE.0) THEN
         JR1 = 0
         JI1 = 0
      ELSE IF (SOPHST) THEN
         JR1 = 14
         JI1 = 4
      ELSE
         JR1 = 5
         JI1 = 1
      END IF
      IF (VECTOL) THEN
         JR2 = NEQF
      ELSE
         JR2 = 1
      END IF
      LRWREQ = 23 + 21*NEQF + JR1*NEQG + JR2*2
      LIWREQ = 21 + JI1*NEQG
      IF (LRWORK.LT.LRWREQ) THEN
         IER = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99982)
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99981) LRWORK, LRWREQ
      ELSE
         IWORK(SVLRWK) = LRWORK
      END IF
      IF (LIWORK.LT.LIWREQ) THEN
         IER = 1
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99980)
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99979) LIWORK, LIWREQ
      ELSE
         IWORK(SVLIWK) = LIWORK
      END IF
      IF (IER.EQ.1) GO TO 100
C
C     CHECK FOR CRITICAL TIME POINT
C
      IF (CRIT) THEN
         IWORK(SVCRIT) = 1
         RWORK(CRITSV) = TCRIT
      ELSE
         IWORK(SVCRIT) = 0
      END IF
C
C     CHECK FOR ONESTEP OR INTERVAL MODE
C
      IF (ONESTP) THEN
         IWORK(SV1STP) = 1
      ELSE
         IWORK(SV1STP) = 0
      END IF
C
C     CHECK FOR MAXIMUM STEP SIZE AND LIMIT ON NUMBER OF STEPS
C
      RWORK(HMAXSV) = ABS(HMAX)
      IF (MAXSTP.GT.0) THEN
         IWORK(SVMXST) = MAXSTP
      ELSE
         IWORK(SVMXST) = 1000
      END IF
C
C     CHECK TOLERANCE VECTORS
C
      DO 20 I = 1, JR2
         IF (ATOL(I).LT.ZERO) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99978) I, ATOL(I)
            GO TO 100
         END IF
   20 CONTINUE
      DO 40 I = 1, JR2
         IF (RTOL(I).LT.ZERO) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99977) I, RTOL(I)
            GO TO 100
         END IF
   40 CONTINUE
      FOURU = FOUR*X02AJF()
      DO 60 I = 1, JR2
         IF (ATOL(I).EQ.ZERO .AND. RTOL(I).LT.FOURU) THEN
            IER = 1
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99976) I
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99975) RTOL(I), FOURU
            GO TO 100
         END IF
   60 CONTINUE
C
C     COPY TOLERANCE VECTORS
C
      IBRTL = 20 + 21*NEQF
      IBATL = JR2 + IBRTL
      DO 80 I = 1, JR2
         RWORK(IBRTL+I) = RTOL(I)
         RWORK(IBATL+I) = ATOL(I)
   80 CONTINUE
C
C     SET IER FOR REFERENCE IN MAIN INTEGRATOR
C
  100 CONTINUE
C
      ERSTAT = IER + 1
      CALL D02QWZ(ERSTAT,SET)
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
C     SET STATEF AND ALTERG TO SUITABLE VALUES FOR A SUBSEQUENT
C     CONTINUATION CALL
C
      IF (IFAIL.EQ.0) THEN
         STATEF = CONTIN
         IF (NEQG.GT.0) ALTERG = .FALSE.
      END IF
C
      RETURN
C
99999 FORMAT (' ** On the first entry STATEF should be ''S'', but is '''
     *       ,A1,'''.')
99998 FORMAT (' ** The value of NEQF,',I16,', is less than 1.')
99997 FORMAT (' ** VECTOL = .TRUE. and the dimension of RTOL is too sm',
     *       'all:')
99996 FORMAT (' ** VECTOL = .TRUE. and the dimension of ATOL is too sm',
     *       'all:')
99995 FORMAT (' ** VECTOL = .FALSE. and the dimension of RTOL is too s',
     *       'mall:')
99994 FORMAT (' ** VECTOL = .FALSE. and the dimension of ATOL is too s',
     *       'mall:')
99993 FORMAT ('    length given (LRTOL) =',I16,', length required is ',
     *       I16,'.')
99992 FORMAT ('    length given (LATOL) =',I16,', length required is ',
     *       I16,'.')
99991 FORMAT ('    length given (LRTOL) =',I16,
     *       ', length required is 1.')
99990 FORMAT ('    length given (LATOL) =',I16,
     *       ', length required is 1.')
99989 FORMAT (' ** On a call with STATEF = ''',A1,''' the value of NEQ',
     *       'F has been changed:')
99988 FORMAT ('    NEQF =',I16,',  and was',I16,'.')
99987 FORMAT (' ** On a call with STATEF = ''',A1,''' an illegal chang',
     *       'e from')
99986 FORMAT ('    scalar tolerences to vector tolerances was specified'
     *       )
99985 FORMAT ('    vector tolerences to scalar tolerances was specified'
     *       )
99984 FORMAT (' ** On a call with STATEF = ''',A1,''' ALTERG = .FALSE.',
     *       ' but NEQG has been')
99983 FORMAT ('    changed to',I16,' from',I16,'.')
99982 FORMAT (' ** With the given input values the dimension of RWORK ',
     *       'is too small:')
99981 FORMAT ('    length given (LRWORK) =',I16,', length required is',
     *       I16,'.')
99980 FORMAT (' ** With the given input values the dimension of IWORK ',
     *       'is too small:')
99979 FORMAT ('    length given (LIWORK) =',I16,', length required is',
     *       I16,'.')
99978 FORMAT (' ** Component',I16,' of ATOL,',1P,D13.5,', is negative.')
99977 FORMAT (' ** Component',I16,' of RTOL,',1P,D13.5,', is negative.')
99976 FORMAT (' ** Component',I16,' of ATOL is 0.0E0 and the correspon',
     *       'ing component')
99975 FORMAT ('    of RTOL,',1P,D13.5,', is less than 4.0*X02AJF(),',1P,
     *       D13.5,'.')
99974 FORMAT (' ** Illegal value for STATEF, ''',A1,'''. Should be ''S',
     *       ''', ''R'' or ''C''.')
      END
