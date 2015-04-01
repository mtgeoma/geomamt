      SUBROUTINE E04NFW(N,NCLIN,TITLE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1595 (JUN 1995).
C
C     ******************************************************************
C     E04NFW loads the default values of parameters not set by the user.
C
C     Original Fortran 77 version written 30-December-1986.
C     This version of E04NFW dated  15-Apr-94.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, HUNDRD
      PARAMETER         (ZERO=0.0D+0,HUNDRD=1.0D+2)
      DOUBLE PRECISION  RDUMMY
      INTEGER           IDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0,IDUMMY=-11111)
      DOUBLE PRECISION  GIGANT
      PARAMETER         (GIGANT=1.0D+20*0.99999D+0)
      DOUBLE PRECISION  WRKTOL
      PARAMETER         (WRKTOL=1.0D-2)
C     .. Scalar Arguments ..
      INTEGER           N, NCLIN
      CHARACTER*(*)     TITLE
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, EPSPT3, EPSPT5,
     *                  EPSPT8, EPSPT9, TOLACT, TOLFEA, TOLINC, TOLRNK,
     *                  TOLX0
      INTEGER           IDEGEN, IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1,
     *                  ITMAX2, ITNFIX, KCHK, KCYCLE, KDEGEN, LCRASH,
     *                  LINES1, LINES2, LPROB, MAXACT, MAXNZ, MINSUM,
     *                  MM, MSGLC, MXFREE, NDEGEN, NN, NNCLIN, NOUT,
     *                  NPROB
      LOGICAL           NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(MXPARM), WMACH(15)
      INTEGER           IPADLC(15), IPSVLC(MXPARM), NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSMCH
      INTEGER           J, LENT, MSGLVL
      CHARACTER*16      KEY
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
      CHARACTER*3       NOYES(0:1)
      CHARACTER*4       ICRSH(0:2)
      CHARACTER*7       QPTYPE(1:10)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          A00AAF, DCOPY, E04NFN, F06DFF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04MF/NEWOPT
      COMMON            /CE04MF/TOLX0, TOLINC, IDEGEN, KDEGEN, NDEGEN,
     *                  ITNFIX, NFIX
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
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
      SAVE              /AX02ZA/, /BE04MF/, /FE04MF/, /GE04MF/
C     .. Data statements ..
      DATA              NOYES(0), NOYES(1)/' NO', 'YES'/
      DATA              ICRSH(0), ICRSH(1), ICRSH(2)/'COLD', 'WARM',
     *                  'HOT '/
      DATA              QPTYPE(1), QPTYPE(2)/'     FP', '     LP'/
      DATA              QPTYPE(3), QPTYPE(4), QPTYPE(5),
     *                  QPTYPE(6)/'    QP1', '    QP2', '    QP3',
     *                  '    QP4'/
      DATA              QPTYPE(7), QPTYPE(8), QPTYPE(9),
     *                  QPTYPE(10)/'       ', '       ', '       ',
     *                  'ILLEGAL'/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
C
C     Make a dummy call to E04NFN to ensure that the defaults are set.
C
      CALL E04NFN(NOUT,'*',KEY)
      NEWOPT = .TRUE.
C
C     Save the optional parameters set by the user.  The values in
C     RPRMLC and IPRMLC may be changed to their default values.
C
      CALL F06DFF(MXPARM,IPRMLC,1,IPSVLC,1)
      CALL DCOPY(MXPARM,RPRMLC,1,RPSVLC,1)
C
      IF (MSGLVL.EQ.IDUMMY) MSGLVL = 10
      IF (IPRNT.LT.0) IPRNT = NOUT
      IF (ISUMRY.LT.0 .OR. MSGLVL.LT.5) ISUMRY = -1
      IPRINT = IPRNT
      ISUMM = ISUMRY
      IF (KCHK.LE.0) KCHK = 50
      IF (KCYCLE.LE.0) KCYCLE = 5
      IF (KCYCLE.GT.9999999) KCYCLE = 9999999
      KDEGEN = KCYCLE
      IF (LPROB.LT.0) LPROB = 4
      IF (LCRASH.LT.0 .OR. LCRASH.GT.2) LCRASH = 0
      IF (ITMAX1.LT.0) ITMAX1 = MAX(50,5*(N+NCLIN))
      IF (ITMAX2.LT.0) ITMAX2 = MAX(50,5*(N+NCLIN))
C
      IF (MAXACT.LT.0 .OR. MAXACT.GT.N .OR. MAXACT.GT.NCLIN)
     *    MAXACT = MAX(1,MIN(N,NCLIN))
C
      IF (MAXNZ.LT.0 .AND. MM.GE.0) MAXNZ = MM
      IF (MAXNZ.EQ.0) MAXNZ = 1
      IF (MAXNZ.LT.1 .OR. MAXNZ.GT.N) MAXNZ = N
      IF (MXFREE.LT.1 .OR. MXFREE.GT.N) MXFREE = N
      IF (MXFREE.LT.MAXNZ) MXFREE = MAXNZ
      IF (MINSUM.LT.0) MINSUM = 0
      IF (NCLIN.LT.N .AND. LPROB.LE.2) THEN
         MXFREE = NCLIN + 1
         MAXNZ = MXFREE
      END IF
      IF (MM.LT.0 .OR. MM.GT.N) MM = N
      IF (LPROB.LE.2) MM = 0
C
      IF (TOLACT.LT.ZERO) TOLACT = WRKTOL
      IF (TOLFEA.EQ.RDUMMY .OR. (TOLFEA.GE.ZERO .AND. TOLFEA.LT.EPSMCH))
     *    TOLFEA = EPSPT5
      IF (TOLRNK.LE.ZERO) TOLRNK = EPSMCH*HUNDRD
      IF (BIGBND.LE.ZERO) BIGBND = GIGANT
      IF (BIGDX.LE.ZERO) BIGDX = MAX(GIGANT,BIGBND)
C
      IF (MSGLVL.GT.0) THEN
C
C        Print the title.
C
         LENT = LEN(TITLE)
         WRITE (REC,FMT=99993) (TITLE(J:J),J=1,LENT)
         CALL X04BAY(IPRINT,2,REC)
         CALL A00AAF
C
         IF (MSGLVL.GE.5 .AND. ISUMM.GE.0 .AND. ISUMM.NE.IPRINT) THEN
            WRITE (REC,FMT=99992) (TITLE(J:J),J=1,LENT)
            CALL X04BAY(ISUMM,2,REC)
         END IF
C
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99998) QPTYPE(LPROB)
         CALL X04BAY(IPRINT,2,REC)
         WRITE (REC,FMT=99997) NCLIN, TOLFEA, N, TOLACT, MM, TOLRNK
         CALL X04BAY(IPRINT,4,REC)
         WRITE (REC,FMT=99996) BIGBND, ICRSH(LCRASH), BIGDX, EPSMCH
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99995) KCHK, KDEGEN, NOYES(MINSUM)
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99994) MAXNZ, MSGLVL, ITMAX1, ISUMRY, ITMAX2
         CALL X04BAY(IPRINT,4,REC)
      END IF
C
      RETURN
C
C
C     End of  E04NFW.  (QPDFLT)
C
99999 FORMAT (/' Parameters',/' ----------')
99998 FORMAT (/' Problem type...........',3X,A7)
99997 FORMAT (/' Linear constraints.....',I10,7X,'Feasibility toleranc',
     *       'e..',1P,D10.2,/' Variables..............',I10,7X,'Crash ',
     *       'tolerance........',1P,D10.2,/' Hessian rows...........',
     *       I10,7X,'Rank tolerance.........',1P,D10.2)
99996 FORMAT (/' Infinite bound size....',1P,D10.2,7X,A4,' start......',
     *       '.......',/' Infinite step size.....',1P,D10.2,7X,'EPS (m',
     *       'achine precision)',1P,D10.2)
99995 FORMAT (/' Check frequency........',I10,7X,'Expand frequency....',
     *       '...',I10,/' Minimum sum of infeas..',7X,A3)
99994 FORMAT (/' Max degrees of freedom.',I10,7X,'Print level.........',
     *       '...',I10,/' Feasibility phase itns.',I10,7X,'Monitoring ',
     *       'file........',I10,/' Optimality  phase itns.',I10)
99993 FORMAT (/80A1)
99992 FORMAT (/11A1,' monitoring information ')
      END
