      SUBROUTINE E04DGS(NOUT,BUFFER,KEY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     E04DGS  decodes the option contained in  BUFFER  in order to set
C     a parameter value in the relevant element of  IPRMCG  or  RPRMCG.
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
C     E04DGS  calls E04UDX and the subprograms
C                 LOOKUP, SCANNR, TOKENS, UPCASE
C     (now called E04UDY, E04UDW, E04UDV, E04UDU)
C     supplied by Informatics General, Inc., Palo Alto, California.
C
C
C     -- Written on 5-June-1986.
C        Sven Hammarling and Janet Welding, NAG Central Office.
C
C     .. Parameters ..
      INTEGER           MXPARM, NIPARM, NRPARM
      PARAMETER         (MXPARM=30,NIPARM=8,NRPARM=5)
      DOUBLE PRECISION  ZERO, RDUMMY
      PARAMETER         (ZERO=0.0D+0,RDUMMY=-11111.0D0)
      INTEGER           MAXKEY, MAXTIE
      PARAMETER         (MAXKEY=21,MAXTIE=9)
      INTEGER           MAXTOK
      PARAMETER         (MAXTOK=5)
      INTEGER           IDUMMY
      PARAMETER         (IDUMMY=-11111)
      LOGICAL           SORTED
      PARAMETER         (SORTED=.TRUE.)
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*16      KEY
      CHARACTER*(*)     BUFFER
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGDX, EPSRF, ETA, FGUESS, FTOL
      INTEGER           IDBGCG, ITMAX, JVRFY1, JVRFY2, LDBGCG, LVERFY,
     *                  MSGCG, NN
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADCG(MXPARM-NRPARM), RPSVCG(MXPARM)
      INTEGER           IPADCG(MXPARM-NIPARM), IPSVCG(MXPARM)
C     .. Local Scalars ..
      DOUBLE PRECISION  RVALUE
      INTEGER           I, IVALUE, LENBUF, LOC1, LOC2, NTOKEN
      LOGICAL           FIRST, MORE, NUMBER
      CHARACTER*16      KEY2, KEY3, VALUE
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMCG(MXPARM)
      INTEGER           IPRMCG(MXPARM)
      CHARACTER*16      KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK)
C     .. External Functions ..
      LOGICAL           E04UDX
      EXTERNAL          E04UDX
C     .. External Subroutines ..
      EXTERNAL          E04UDV, E04UDY, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         INDEX, LEN
C     .. Common blocks ..
      COMMON            /AE04DG/IPSVCG, IDBGCG, ITMAX, JVRFY1, JVRFY2,
     *                  LDBGCG, LVERFY, MSGCG, NN, IPADCG
      COMMON            /BE04DG/RPSVCG, BIGDX, EPSRF, ETA, FGUESS, FTOL,
     *                  RPADCG
C     .. Equivalences ..
      EQUIVALENCE       (IPRMCG(1),IDBGCG), (RPRMCG(1),BIGDX)
C     .. Save statement ..
      SAVE              /AE04DG/, /BE04DG/, FIRST
C     .. Data statements ..
      DATA              FIRST/.TRUE./
      DATA              KEYS/'BEGIN           ', 'DEBUG           ',
     *                  'DEFAULTS        ', 'END             ',
     *                  'ESTIMATED       ', 'FUNCTION        ',
     *                  'INFINITE        ', 'IPRMCG          ',
     *                  'ITERATIONS      ', 'ITERS:ITERATIONS',
     *                  'ITNS :ITERATIONS', 'LINESEARCH      ',
     *                  'LIST            ', 'MAXIMUM         ',
     *                  'NOLIST          ', 'OPTIMALITY      ',
     *                  'PRINT           ', 'RPRMCG          ',
     *                  'START           ', 'STOP            ',
     *                  'VERIFY          '/
      DATA              TIES/'GRADIENTS       ', 'LEVEL           ',
     *                  'NO              ', 'NO.      :NUMBER',
     *                  'NUMBER          ', 'OBJECTIVE       ',
     *                  'START           ', 'TOLERANCE       ',
     *                  'YES             '/
C     .. Executable Statements ..
C
C     First time in initialize the optional parameters with dummy
C     values.
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         DO 20 I = 1, MXPARM
            IPRMCG(I) = IDUMMY
            RPRMCG(I) = RDUMMY
   20    CONTINUE
      END IF
C
C     Eliminate comments and empty lines.
C     A '*' (the comment symbol) appearing anywhere in BUFFER terminates
C     the string.
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
C     Extract up to MAXTOK tokens from the record.
C     NTOKEN returns how many were actually found.
C     KEY, KEY2, KEY3 are the first tokens if any, otherwise blank.
C
      CALL E04UDV(BUFFER(1:LENBUF),MAXTOK,NTOKEN,TOKEN)
      KEY = TOKEN(1)
      KEY2 = TOKEN(2)
      KEY3 = TOKEN(3)
C
C     Certain keywords require no action.
C
      IF ((KEY.EQ.' ') .OR. (KEY.EQ.'BEGIN') .OR. (KEY.EQ.'LIST')
     *    .OR. (KEY.EQ.'NOLIST') .OR. (KEY.EQ.'END')) GO TO 80
C
C     Most keywords will have an associated integer or real value,
C     so look for it no matter what the keyword.
C
      I = 1
      NUMBER = .FALSE.
C
   40 IF ((I.LT.NTOKEN) .AND. ( .NOT. NUMBER)) THEN
         I = I + 1
         VALUE = TOKEN(I)
         NUMBER = E04UDX(VALUE)
         GO TO 40
      END IF
C
      IF (NUMBER) THEN
         READ (VALUE,FMT='( BN, E16.0 )') RVALUE
      ELSE
         RVALUE = ZERO
      END IF
C
C     Convert the keywords to their most fundamental form
C     (upper case, no abbreviations).
C     SORTED says whether the dictionaries are in alphabetic order.
C     LOCi   says where the keywords are in the dictionaries.
C     LOCi = 0 signals that the keyword wasn't there.
C
      CALL E04UDY(MAXKEY,KEYS,SORTED,KEY,LOC1)
      CALL E04UDY(MAXTIE,TIES,SORTED,KEY2,LOC2)
C
C     Decide what to do about each keyword.
C     The second keyword (if any) might be needed to break ties.
C     Some seemingly redundant testing of MORE is used
C     to avoid compiler limits on the number of consecutive ELSE IFs.
C
      MORE = .TRUE.
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'DEBUG       ') THEN
            IF (KEY2.EQ.'LEVEL       ') THEN
               LDBGCG = RVALUE
            ELSE IF (KEY2.EQ.'START       ') THEN
               IDBGCG = RVALUE
            ELSE IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'DEFAULTS    ') THEN
            DO 60 I = 1, MXPARM
               IPRMCG(I) = IDUMMY
               RPRMCG(I) = RDUMMY
   60       CONTINUE
         ELSE IF (KEY.EQ.'ESTIMATED   ') THEN
            FGUESS = RVALUE
         ELSE IF (KEY.EQ.'FUNCTION    ') THEN
            EPSRF = RVALUE
         ELSE
            MORE = .TRUE.
         END IF
      END IF
C
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'IPRMCG      ') THEN
C           Allow things like  IPRMCG 21 = 37  to set  IPRMCG( 21 ) = 37
            IVALUE = RVALUE
            IF ((IVALUE.GE.1) .AND. (IVALUE.LE.MXPARM)) THEN
               READ (KEY3,FMT='( BN, I16 )') IPRMCG(IVALUE)
            ELSE
               WRITE (REC,FMT=99997) IVALUE
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'ITERATIONS  ') THEN
            ITMAX = RVALUE
         ELSE IF (KEY.EQ.'LINESEARCH  ') THEN
            ETA = RVALUE
         ELSE IF (KEY.EQ.'MAXIMUM     ') THEN
            BIGDX = RVALUE
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
            MSGCG = RVALUE
         ELSE IF (KEY.EQ.'RPRMCG      ') THEN
C           Allow things like  RPRMCG 21 = 2  to set  RPRMCG( 21 ) = 2.0
            IVALUE = RVALUE
            IF ((IVALUE.GE.1) .AND. (IVALUE.LE.MXPARM)) THEN
               READ (KEY3,FMT='( BN, E16.0 )') RPRMCG(IVALUE)
            ELSE
               WRITE (REC,FMT=99997) IVALUE
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'START       ') THEN
            JVRFY1 = RVALUE
         ELSE IF (KEY.EQ.'STOP        ') THEN
            JVRFY2 = RVALUE
         ELSE IF (KEY.EQ.'VERIFY      ') THEN
            IF (KEY2.EQ.'GRADIENTS   ') THEN
               LVERFY = 1
            ELSE IF (KEY2.EQ.'LEVEL       ') THEN
               LVERFY = RVALUE
            ELSE IF (KEY2.EQ.'NO          ') THEN
               LVERFY = -1
            ELSE IF (KEY2.EQ.'OBJECTIVE   ') THEN
               LVERFY = 1
            ELSE IF (KEY2.EQ.'YES         ') THEN
               LVERFY = 1
            ELSE IF (LOC2.EQ.0) THEN
               LVERFY = 1
            END IF
         ELSE
            WRITE (REC,FMT=99999) KEY
            CALL X04BAF(NOUT,REC)
         END IF
      END IF
C
   80 RETURN
C
C
C     End of  E04DGS. ( CGKEY )
C
99999 FORMAT (' XXX  Keyword not recognized:            ',A)
99998 FORMAT (' XXX  Second keyword not recognized:     ',A)
99997 FORMAT (' XXX  The PARM subscript is out of range:',I10)
      END
