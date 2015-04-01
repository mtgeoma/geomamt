      SUBROUTINE E04UCQ(NOUT,BUFFER,KEY)
C     MARK 13 RE-ISSUE.  NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1085 (JUL 1993).
C     MARK 17 REVISED. IER-1605 (JUN 1995).
C
C     ******************************************************************
C     E04UCQ   decodes the option contained in  BUFFER  in order to set
C     a parameter value in the relevant element of the parameter arrays.
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
C     E04UCQ  calls E04UDX and the subprograms
C                 LOOKUP, SCANNR, TOKENS, UPCASE
C     (now called E04UDY, E04UDW, E04UDV, E04UDU)
C     supplied by Informatics General, Inc., Palo Alto, California.
C
C     Systems Optimization Laboratory, Stanford University.
C     This version of E04UCQ  dated 19-Oct-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      INTEGER           MAXKEY, MAXTIE, MAXTOK
      PARAMETER         (MAXKEY=40,MAXTIE=19,MAXTOK=10)
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
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT, CTOL,
     *                  DXLIM, EPSRF, ETA, FDINT, FTOL, HCNDBD, TOLACT,
     *                  TOLFEA, TOLRNK
      INTEGER           IPRNT, ISUMRY, ITMAX1, ITMAX2, ITMXNP, JVRFY1,
     *                  JVRFY2, JVRFY3, JVRFY4, KSAVE, LCRASH, LFORMH,
     *                  LPROB, LVERFY, LVLDER, MSGLS, MSGNP, NLNF, NLNJ,
     *                  NLNX, NLOAD, NN, NNCLIN, NNCNLN, NPROB, NSAVE
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNP(22), RPSVLS(MXPARM),
     *                  RPSVNP(MXPARM)
      INTEGER           IPADLS(19), IPADNP(15), IPSVLS(MXPARM),
     *                  IPSVNP(MXPARM)
C     .. Local Scalars ..
      DOUBLE PRECISION  RVALUE
      INTEGER           I, IVALUE, J, LENBUF, LOC1, LOC2, MSGQP, NMAJOR,
     *                  NMINOR, NTOKEN
      LOGICAL           FIRST, MORE, NUMBER
      CHARACTER*16      KEY2, KEY3, VALUE
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNP(MXPARM)
      CHARACTER*16      KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK)
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
      COMMON            /GE04UC/IPSVNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, LVLDER, LVERFY, MSGNP, NLNF, NLNJ, NLNX,
     *                  NNCNLN, NSAVE, NLOAD, KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),ITMXNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (ITMXNP,NMAJOR), (ITMAX2,NMINOR), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /DE04NC/, /EE04NC/, /GE04UC/, /HE04UC/, FIRST
C     .. Data statements ..
      DATA              (KEYS(J),J=1,19)/'BEGIN           ',
     *                  'CENTRAL         ', 'COLD            ',
     *                  'CONDITION       ', 'CONSTRAINTS     ',
     *                  'CRASH           ', 'DEFAULTS        ',
     *                  'DERIVATIVE      ', 'DIFFERENCE      ',
     *                  'END             ', 'FEASIBILITY     ',
     *                  'FUNCTION        ', 'HESSIAN         ',
     *                  'INFINITE        ', 'IPRMLS          ',
     *                  'ITERATIONS      ', 'ITERS:ITERATIONS',
     *                  'ITNS :ITERATIONS', 'LINE            '/
      DATA              (KEYS(J),J=20,40)/'LINEAR          ',
     *                  'LINESEARCH:LINE ', 'LIST            ',
     *                  'LOWER           ', 'MAJOR           ',
     *                  'MINOR           ', 'MONITORING      ',
     *                  'NOLIST          ', 'NONLINEAR       ',
     *                  'OPTIMALITY      ', 'PRINT           ',
     *                  'PROBLEM         ', 'ROW             ',
     *                  'RPRMLS          ', 'START           ',
     *                  'STEP            ', 'STOP            ',
     *                  'UPPER           ', 'VARIABLES       ',
     *                  'VERIFY          ', 'WARM            '/
      DATA              FIRST/.TRUE./
      DATA              TIES/'BOUND           ', 'CONSTRAINTS     ',
     *                  'FEASIBILITY     ', 'GRADIENTS       ',
     *                  'ITERATIONS      ', 'ITERS:ITERATIONS',
     *                  'ITNS :ITERATIONS', 'JACOBIAN        ',
     *                  'LEVEL           ', 'NO              ',
     *                  'NO.      :NUMBER', 'NUMBER          ',
     *                  'OBJECTIVE       ', 'PRINT           ',
     *                  'SEARCH          ', 'STEP            ',
     *                  'TOLERANCE       ', 'VARIABLES       ',
     *                  'YES             '/
C     .. Executable Statements ..
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         DO 20 I = 1, MXPARM
            RPRMLS(I) = RDUMMY
            IPRMLS(I) = IDUMMY
            RPRMNP(I) = RDUMMY
            IPRMNP(I) = IDUMMY
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
         IF (KEY.EQ.'CENTRAL     ') THEN
            CDINT = RVALUE
         ELSE IF (KEY.EQ.'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY.EQ.'CONDITION   ') THEN
            HCNDBD = RVALUE
         ELSE IF (KEY.EQ.'CONSTRAINTS ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY.EQ.'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY.EQ.'DEFAULTS    ') THEN
            DO 60 I = 1, MXPARM
               IPRMLS(I) = IDUMMY
               RPRMLS(I) = RDUMMY
               IPRMNP(I) = IDUMMY
               RPRMNP(I) = RDUMMY
   60       CONTINUE
         ELSE IF (KEY.EQ.'DERIVATIVE  ') THEN
            LVLDER = RVALUE
         ELSE IF (KEY.EQ.'DIFFERENCE  ') THEN
            FDINT = RVALUE
         ELSE IF (KEY.EQ.'FEASIBILITY ') THEN
            TOLFEA = RVALUE
            CTOL = RVALUE
         ELSE IF (KEY.EQ.'FUNCTION    ') THEN
            EPSRF = RVALUE
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
            END IF
         ELSE IF (KEY.EQ.'IPRMLS      ') THEN
C           Allow things like  IPRMLS 21 = 100  to set IPRMLS(21) = 100
            IVALUE = RVALUE
            IF (IVALUE.GE.1 .AND. IVALUE.LE.MXPARM) THEN
               READ (KEY3,FMT='(BN, I16)') IPRMLS(IVALUE)
            ELSE
               WRITE (REC,FMT=99997) IVALUE
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'ITERATIONS  ') THEN
            NMAJOR = RVALUE
         ELSE IF (KEY.EQ.'LINE        ') THEN
            ETA = RVALUE
         ELSE IF (KEY.EQ.'LINEAR      ') THEN
            IF (KEY2.EQ.'CONSTRAINTS ') NNCLIN = RVALUE
            IF (KEY2.EQ.'FEASIBILITY ') TOLFEA = RVALUE
            IF (KEY2.EQ.'SEARCH      ') ETA = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'LOWER       ') THEN
            BNDLOW = RVALUE
         ELSE
            MORE = .TRUE.
         END IF
      END IF
C
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'MAJOR       ') THEN
            IF (KEY2.EQ.'ITERATIONS  ') NMAJOR = RVALUE
            IF (KEY2.EQ.'PRINT       ') MSGNP = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'MINOR       ') THEN
            IF (KEY2.EQ.'ITERATIONS  ') NMINOR = RVALUE
            IF (KEY2.EQ.'PRINT       ') MSGQP = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'MONITORING  ') THEN
            ISUMRY = RVALUE
         ELSE IF (KEY.EQ.'NONLINEAR   ') THEN
            IF (KEY2.EQ.'CONSTRAINTS ') NNCNLN = RVALUE
            IF (KEY2.EQ.'FEASIBILITY ') CTOL = RVALUE
            IF (KEY2.EQ.'JACOBIAN    ') NLNJ = RVALUE
            IF (KEY2.EQ.'OBJECTIVE   ') NLNF = RVALUE
            IF (KEY2.EQ.'VARIABLES   ') NLNX = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'OPTIMALITY  ') THEN
            FTOL = RVALUE
         ELSE
            MORE = .TRUE.
         END IF
      END IF
C
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'PRINT       ') THEN
            MSGNP = RVALUE
         ELSE IF (KEY.EQ.'PROBLEM     ') THEN
            IF (KEY2.EQ.'NUMBER      ') NPROB = RVALUE
         ELSE IF (KEY.EQ.'ROW         ') THEN
            IF (KEY2.EQ.'TOLERANCE   ') CTOL = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'RPRMLS      ') THEN
C           Allow things like  RPRMLS 21 = 2  to set RPRMLS(21) = 2.0
            IVALUE = RVALUE
            IF (IVALUE.GE.1 .AND. IVALUE.LE.MXPARM) THEN
               READ (KEY3,FMT='(BN, E16.0)') RPRMLS(IVALUE)
            ELSE
               WRITE (REC,FMT=99997) IVALUE
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'START       ') THEN
            IF (KEY2.EQ.'CONSTRAINTS ') JVRFY3 = RVALUE
            IF (KEY2.EQ.'OBJECTIVE   ') JVRFY1 = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'STEP        ') THEN
            DXLIM = RVALUE
         ELSE IF (KEY.EQ.'STOP        ') THEN
            IF (KEY2.EQ.'CONSTRAINTS ') JVRFY4 = RVALUE
            IF (KEY2.EQ.'OBJECTIVE   ') JVRFY2 = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'UPPER       ') THEN
            BNDUPP = RVALUE
         ELSE IF (KEY.EQ.'VARIABLES   ') THEN
            NN = RVALUE
         ELSE IF (KEY.EQ.'VERIFY      ') THEN
            IF (KEY2.EQ.'OBJECTIVE   ') LVERFY = 1
            IF (KEY2.EQ.'CONSTRAINTS ') LVERFY = 2
            IF (KEY2.EQ.'NO          ') LVERFY = -1
            IF (KEY2.EQ.'YES         ') LVERFY = 3
            IF (KEY2.EQ.'GRADIENTS   ') LVERFY = 3
            IF (KEY2.EQ.'LEVEL       ') LVERFY = RVALUE
            IF (LOC2.EQ.0) LVERFY = 3
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
C     End of  E04UCQ. (NPKEY)
C
99999 FORMAT (' XXX  Keyword not recognized:         ',A)
99998 FORMAT (' XXX  Second keyword not recognized:  ',A)
99997 FORMAT (' XXX  The PARM subscript is out of range:',I10)
99996 FORMAT (' XXX  Ambiguous keyword:              ',A)
      END
