      SUBROUTINE E04MFX(NOUT,BUFFER,KEY)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1568 (JUN 1995).
C
C     ******************************************************************
C     E04MFX   decodes the option contained in  BUFFER  in order to set
C     a parameter value in the relevant element of  IPRMLC  or  RPRMLC.
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
C     E04MFX  calls E04UDX and the subprograms
C                 LOOKUP, SCANNR, TOKENS, UPCASE
C     (now called E04UDY, E04UDW, E04UDV, E04UDU)
C     supplied by Informatics General, Inc., Palo Alto, California.
C
C     This version of E04MFX dated 10-Apr-94.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      INTEGER           MAXKEY, MAXTIE, MAXTOK, MAXTYP
      PARAMETER         (MAXKEY=20,MAXTIE=5,MAXTOK=10,MAXTYP=4)
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
      INTEGER           IPRNT, ISUMRY, ITMAX1, ITMAX2, KCHK, KCYCLE,
     *                  LCRASH, LPROB, MAXACT, MAXNZ, MINSUM, MM, MSGLC,
     *                  MXFREE, NN, NNCLIN, NPROB
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(MXPARM)
      INTEGER           IPADLC(15), IPSVLC(MXPARM)
C     .. Local Scalars ..
      DOUBLE PRECISION  RVALUE
      INTEGER           I, LENBUF, LOC1, LOC2, LOC3, MSGLVL, NTOKEN
      LOGICAL           FIRST, MORE, NUMBER
      CHARACTER*16      KEY2, KEY3, VALUE
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
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
      COMMON            /FE04MF/IPSVLC, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  KCHK, KCYCLE, LCRASH, LPROB, MAXACT, MXFREE,
     *                  MAXNZ, MM, MINSUM, MSGLC, NN, NNCLIN, NPROB,
     *                  IPADLC
      COMMON            /GE04MF/RPSVLC, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLC
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLC(1),IPRNT), (RPRMLC(1),BIGBND)
      EQUIVALENCE       (MSGLC,MSGLVL)
C     .. Save statement ..
      SAVE              /FE04MF/, /GE04MF/, FIRST
C     .. Data statements ..
      DATA              FIRST/.TRUE./
      DATA              KEYS(1)/'BEGIN           '/,
     *                  KEYS(2)/'CHECK           '/,
     *                  KEYS(3)/'COLD            '/,
     *                  KEYS(4)/'CRASH           '/,
     *                  KEYS(5)/'DEFAULTS        '/,
     *                  KEYS(6)/'END             '/,
     *                  KEYS(7)/'EXPAND          '/,
     *                  KEYS(8)/'FEASIBILITY     '/,
     *                  KEYS(9)/'INFINITE        '/,
     *                  KEYS(10)/'ITERATIONS      '/,
     *                  KEYS(11)/'ITERS:ITERATIONS'/,
     *                  KEYS(12)/'ITNS :ITERATIONS'/,
     *                  KEYS(13)/'LIST            '/,
     *                  KEYS(14)/'MIN     :MINIMIZ'/,
     *                  KEYS(15)/'MINIMIZ         '/,
     *                  KEYS(16)/'MONITORING      '/,
     *                  KEYS(17)/'NOLIST          '/,
     *                  KEYS(18)/'PRINT           '/,
     *                  KEYS(19)/'PROBLEM         '/,
     *                  KEYS(20)/'WARM            '/
      DATA              TIES/'BOUND           ', 'NO              ',
     *                  'STEP            ', 'TYPE            ',
     *                  'YES             '/
      DATA              TYPE/'FEASIBLE     :FP', 'FP              ',
     *                  'LINEAR       :LP', 'LP              '/
C     .. Executable Statements ..
C     ------------------------------------------------------------------
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         DO 20 I = 1, MXPARM
            IPRMLC(I) = IDUMMY
            RPRMLC(I) = RDUMMY
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
      NTOKEN = MAXTOK
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
         IF (KEY.EQ.'CHECK       ') THEN
            KCHK = RVALUE
         ELSE IF (KEY.EQ.'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY.EQ.'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY.EQ.'DEFAULTS    ') THEN
            DO 60 I = 1, MXPARM
               IPRMLC(I) = IDUMMY
               RPRMLC(I) = RDUMMY
   60       CONTINUE
         ELSE IF (KEY.EQ.'EXPAND      ') THEN
            KCYCLE = RVALUE
         ELSE IF (KEY.EQ.'FEASIBILITY ') THEN
            TOLFEA = RVALUE
         ELSE
            MORE = .TRUE.
         END IF
      END IF
C
      IF (MORE) THEN
         MORE = .FALSE.
         IF (KEY.EQ.'INFINITE    ') THEN
            IF (KEY2.EQ.'BOUND       ') BIGBND = RVALUE*0.99999D+0
            IF (KEY2.EQ.'STEP        ') BIGDX = RVALUE
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'ITERATIONS  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY.EQ.'MINIMIZ     ') THEN
            IF (KEY3.EQ.'YES         ') MINSUM = 1
            IF (KEY3.EQ.'NO          ') MINSUM = 0
            IF (LOC2.EQ.0) THEN
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
         ELSE IF (KEY.EQ.'MONITORING  ') THEN
            ISUMRY = RVALUE
         ELSE IF (KEY.EQ.'PRINT       ') THEN
            MSGLVL = RVALUE
         ELSE IF (KEY.EQ.'PROBLEM     ') THEN
            IF (KEY2.EQ.'TYPE  ') THEN
C
C              Recognize     Problem type = LP     etc.
C
               CALL E04UDY(MAXTYP,TYPE,SORTED,KEY3,LOC3)
               IF (KEY3.EQ.'FP') LPROB = 1
               IF (KEY3.EQ.'LP') LPROB = 2
               IF (LOC3.EQ.0) THEN
                  WRITE (REC,FMT=99997) KEY3
                  CALL X04BAF(NOUT,REC)
                  LPROB = 10
               END IF
            ELSE
               WRITE (REC,FMT=99998) KEY2
               CALL X04BAF(NOUT,REC)
            END IF
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
C     End of  E04MFX.  (LPKEY)
C
99999 FORMAT (' XXX  Keyword not recognized:         ',A)
99998 FORMAT (' XXX  Second keyword not recognized:  ',A)
99997 FORMAT (' XXX  Third  keyword not recognized:  ',A)
99996 FORMAT (' XXX  Ambiguous keyword:              ',A)
      END
