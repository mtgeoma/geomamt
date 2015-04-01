      SUBROUTINE E04NCN(NOUT,BUFFER,KEY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-1066 (JUL 1993).
C     MARK 17 REVISED. IER-1578 (JUN 1995).
C
C     ******************************************************************
C     E04NCN   decodes the option contained in  BUFFER  in order to set
C     a parameter value in the relevant element of  IPRMLS  or  RPRMLS.
C
C
C     Input:
C
C     NOUT   A unit number for printing error messages.
C            NOUT  must be a valid unit.
C
C     Output:
C
C     KEY    The first keyword contained in BUFFER.
C
C
C     E04NCN  calls E04UDX and the subprograms
C                 LOOKUP, SCANNR, TOKENS, UPCASE
C     (now called E04UDY, E04UDW, E04UDV, E04UDU)
C     supplied by Informatics General, Inc., Palo Alto, California.
C
C     Systems Optimization Laboratory, Stanford University.
C     This version dated 21-Apr-93.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      INTEGER           MAXKEY, MAXTIE, MAXTOK, MAXTYP
      PARAMETER         (MAXKEY=24,MAXTIE=12,MAXTOK=10,MAXTYP=16)
      INTEGER           IDUMMY
      DOUBLE PRECISION  RDUMMY
      LOGICAL           SORTED
      DOUBLE PRECISION  ZERO
      PARAMETER         (IDUMMY=-11111,RDUMMY=-11111.0D+0,SORTED=.TRUE.,
     *                  ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*16      KEY
      CHARACTER*(*)     BUFFER
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, TOLACT, TOLFEA,
     *                  TOLRNK
      INTEGER           IPRNT, ISUMRY, ITMAX1, ITMAX2, LCRASH, LFORMH,
     *                  LPROB, MSGLS, NN, NNCLIN, NPROB
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPSVLS(MXPARM)
      INTEGER           IPADLS(19), IPSVLS(MXPARM)
C     .. Local Scalars ..
      DOUBLE PRECISION  RVALUE
      INTEGER           I, LENBUF, LOC1, LOC2, LOC3, NTOKEN
      LOGICAL           FIRST, MORE, NUMBER
      CHARACTER*16      KEY2, KEY3, VALUE
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM)
      INTEGER           IPRMLS(MXPARM)
      CHARACTER*16      KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK),
     *                  TYPE(MAXTYP)
C     .. External Functions ..
      LOGICAL           E04UDX
      EXTERNAL          E04UDX
C     .. External Subroutines ..
      EXTERNAL          E04UDV, E04UDY, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         INDEX, LEN
C     .. Common blocks ..
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
C     .. Save statement ..
      SAVE              /DE04NC/, /EE04NC/, FIRST
C     .. Data statements ..
      DATA              FIRST/.TRUE./
      DATA              KEYS/'BEGIN           ', 'COLD            ',
     *                  'CONSTRAINTS     ', 'CRASH           ',
     *                  'DEFAULTS        ', 'END             ',
     *                  'FEASIBILITY     ', 'HESSIAN         ',
     *                  'INFINITE        ', 'ITERATIONS      ',
     *                  'ITERS:ITERATIONS', 'ITNS :ITERATIONS',
     *                  'LINEAR          ', 'LIST            ',
     *                  'LOWER           ', 'MONITORING      ',
     *                  'NOLIST          ', 'OPTIMALITY      ',
     *                  'PRINT           ', 'PROBLEM         ',
     *                  'RANK            ', 'UPPER           ',
     *                  'VARIABLES       ', 'WARM            '/
      DATA              TIES/'BOUND           ', 'CONSTRAINTS     ',
     *                  'FILE            ', 'LEVEL           ',
     *                  'NO              ', 'NO.      :NUMBER',
     *                  'NUMBER          ', 'PHASE           ',
     *                  'STEP            ', 'TOLERANCE       ',
     *                  'TYPE            ', 'YES             '/
      DATA              TYPE/'FP              ', 'LEAST       :LS1',
     *                  'LINEAR       :LP', 'LP              ',
     *                  'LS          :LS1', 'LS1             ',
     *                  'LS2             ', 'LS3             ',
     *                  'LS4             ', 'LSQ         :LS1',
     *                  'QP          :QP2', 'QP1             ',
     *                  'QP2             ', 'QP3             ',
     *                  'QP4             ', 'QUADRATIC   :QP2'/
C     .. Executable Statements ..
      IF (FIRST) THEN
         FIRST = .FALSE.
         DO 20 I = 1, MXPARM
            IPRMLS(I) = IDUMMY
            RPRMLS(I) = RDUMMY
   20    CONTINUE
      END IF
C
C     Eliminate comments and empty lines.
C     A '*' appearing anywhere in BUFFER terminates the string.
C
      I = INDEX(BUFFER,'*')
      IF (I.EQ.0) THEN
         LENBUF = LEN(BUFFER)
      ELSE
         LENBUF = I - 1
      END IF
      IF (LENBUF.LE.0) THEN
         KEY = '*'
         GO TO 80
      END IF
C
C     ------------------------------------------------------------------
C     Extract up to MAXTOK tokens from the record.
C     NTOKEN returns how many were actually found.
C     KEY, KEY2, KEY3 are the first tokens if any, otherwise blank.
C     ------------------------------------------------------------------
      CALL E04UDV(BUFFER(1:LENBUF),MAXTOK,NTOKEN,TOKEN)
      KEY = TOKEN(1)
      KEY2 = TOKEN(2)
      KEY3 = TOKEN(3)
C
C     Certain keywords require no action.
C
      IF (KEY.EQ.' ' .OR. KEY.EQ.'BEGIN') GO TO 80
      IF (KEY.EQ.'LIST' .OR. KEY.EQ.'NOLIST') GO TO 80
      IF (KEY.EQ.'END') GO TO 80
C
C     Most keywords will have an associated integer or real value,
C     so look for it no matter what the keyword.
C
      I = 1
      NUMBER = .FALSE.
C
   40 IF (I.LT.NTOKEN .AND. .NOT. NUMBER) THEN
         I = I + 1
         VALUE = TOKEN(I)
         NUMBER = E04UDX(VALUE)
         GO TO 40
      END IF
C
      IF (NUMBER) THEN
         READ (VALUE,FMT='(BN, E16.0)') RVALUE
      ELSE
         RVALUE = ZERO
      END IF
C
C     Convert the keywords to their most fundamental form
C     (upper case, no abbreviations).
C     SORTED says whether the dictionaries are in alphabetic order.
C     LOCi   says where the keywords are in the dictionaries.
C     LOCi = 0 signals that the keyword wasn't there.
C     LOCi < 0 signals that the keyword is ambiguous.
C
      CALL E04UDY(MAXKEY,KEYS,SORTED,KEY,LOC1)
      IF (LOC1.LT.0) THEN
         WRITE (REC,FMT=99996) KEY
         CALL X04BAF(NOUT,REC)
         RETURN
      END IF
      CALL E04UDY(MAXTIE,TIES,SORTED,KEY2,LOC2)
C
C     ------------------------------------------------------------------
C     Decide what to do about each keyword.
C     The second keyword (if any) might be needed to break ties.
C     Some seemingly redundant testing of MORE is used
C     to avoid compiler limits on the number of consecutive ELSE IFs.
C     ------------------------------------------------------------------
      MORE = .TRUE.
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY.EQ.'CONSTRAINTS ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY.EQ.'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY.EQ.'DEFAULTS    ') THEN
            DO 60 I = 1, MXPARM
               IPRMLS(I) = IDUMMY
               RPRMLS(I) = RDUMMY
   60       CONTINUE
         ELSE IF (KEY.EQ.'FEASIBILITY ') THEN
            IF (KEY2.EQ.'PHASE       ') ITMAX1 = RVALUE
            IF (KEY2.EQ.'TOLERANCE   ') TOLFEA = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            ELSE IF (LOC2.LT.0) THEN
               WRITE (REC,FMT=99996) KEY2
               CALL X04BAF(NOUT,REC)
               RETURN
            END IF
         ELSE
            MORE = .TRUE.
         END IF
      END IF
C
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'HESSIAN     ') THEN
            LFORMH = 0
            IF (KEY2.EQ.'YES         ') LFORMH = 1
         ELSE IF (KEY.EQ.'INFINITE    ') THEN
            IF (KEY2.EQ.'BOUND       ') BIGBND = RVALUE*0.99999D+0
            IF (KEY2.EQ.'STEP        ') BIGDX = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            ELSE IF (LOC2.LT.0) THEN
               WRITE (REC,FMT=99996) KEY2
               CALL X04BAF(NOUT,REC)
               RETURN
            END IF
         ELSE IF (KEY.EQ.'ITERATIONS  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY.EQ.'LINEAR      ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY.EQ.'LOWER       ') THEN
            BNDLOW = RVALUE
         ELSE
            MORE = .TRUE.
         END IF
      END IF
C
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'MONITORING  ') THEN
            ISUMRY = RVALUE
         ELSE IF (KEY.EQ.'OPTIMALITY  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY.EQ.'PROBLEM     ') THEN
            IF (KEY2.EQ.'NUMBER') THEN
               NPROB = RVALUE
            ELSE IF (KEY2.EQ.'TYPE  ') THEN
C
C              Recognize     Problem type = LP     etc.
C
               CALL E04UDY(MAXTYP,TYPE,SORTED,KEY3,LOC3)
               IF (KEY3.EQ.'FP') LPROB = 1
               IF (KEY3.EQ.'LP') LPROB = 2
               IF (KEY3.EQ.'QP1') LPROB = 3
               IF (KEY3.EQ.'QP2') LPROB = 4
               IF (KEY3.EQ.'QP3') LPROB = 5
               IF (KEY3.EQ.'QP4') LPROB = 6
               IF (KEY3.EQ.'LS1') LPROB = 7
               IF (KEY3.EQ.'LS2') LPROB = 8
               IF (KEY3.EQ.'LS3') LPROB = 9
               IF (KEY3.EQ.'LS4') LPROB = 10
               IF (LOC3.EQ.0) THEN
                  WRITE (REC,FMT=99997) KEY3
                  CALL X04BAF(NOUT,REC)
               END IF
            ELSE
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE
            MORE = .TRUE.
         END IF
      END IF
C
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'PRINT       ') THEN
            IF (KEY2.EQ.'FILE        ') IPRNT = RVALUE
            IF (KEY2.EQ.'LEVEL       ') MSGLS = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998)
               CALL X04BAF(NOUT,REC)
            ELSE IF (LOC2.LT.0) THEN
               WRITE (REC,FMT=99996) KEY2
               CALL X04BAF(NOUT,REC)
               RETURN
            END IF
         ELSE IF (KEY.EQ.'RANK        ') THEN
            TOLRNK = RVALUE
         ELSE IF (KEY.EQ.'UPPER       ') THEN
            BNDUPP = RVALUE
         ELSE IF (KEY.EQ.'VARIABLES   ') THEN
            NN = RVALUE
         ELSE IF (KEY.EQ.'WARM        ') THEN
            LCRASH = 1
         ELSE
            WRITE (REC,FMT=99999) KEY
            CALL X04BAF(NOUT,REC)
         END IF
      END IF
C
   80 RETURN
C
C
C     End of  E04NCN. (LSKEY)
C
99999 FORMAT (' XXX  Keyword not recognized:         ',A)
99998 FORMAT (' XXX  Second keyword not recognized:  ',A)
99997 FORMAT (' XXX  Third  keyword not recognized:  ',A)
99996 FORMAT (' XXX  Ambiguous keyword:              ',A)
      END
