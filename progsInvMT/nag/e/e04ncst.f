      SUBROUTINE E04NCS(M,N,NCLIN,TITLE)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-598 (MAR 1988).
C     MARK 13 REVISED. IER-644 (APR 1988).
C     MARK 16 REVISED. IER-1070 (JUL 1993).
C     MARK 17 REVISED. IER-1582 (JUN 1995).
C
C     ******************************************************************
C     E04NCS  loads the default values of parameters not set by the
C     user.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written 17-September-1985.
C     This version of E04NCS dated  22-Mar-93.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, TEN
      PARAMETER         (ZERO=0.0D+0,TEN=10.0D+0)
      DOUBLE PRECISION  HUNDRD
      PARAMETER         (HUNDRD=100.0D+0)
      DOUBLE PRECISION  RDUMMY
      INTEGER           IDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0,IDUMMY=-11111)
      DOUBLE PRECISION  GIGANT
      PARAMETER         (GIGANT=1.0D+20*0.99999D+0)
      DOUBLE PRECISION  WRKTOL
      PARAMETER         (WRKTOL=1.0D-2)
C     .. Scalar Arguments ..
      INTEGER           M, N, NCLIN
      CHARACTER*(*)     TITLE
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, EPSPT3, EPSPT5,
     *                  EPSPT8, EPSPT9, TOLACT, TOLFEA, TOLRNK
      INTEGER           IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LINES1, LINES2, LPROB, MSGLS,
     *                  NN, NNCLIN, NOUT, NPROB
      LOGICAL           NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPSVLS(MXPARM), WMACH(15)
      INTEGER           IPADLS(19), IPSVLS(MXPARM)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSMCH
      INTEGER           J, LENT, MSGLVL
      CHARACTER*16      KEY
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM)
      INTEGER           IPRMLS(MXPARM)
      CHARACTER*3       CHESS(0:1), LSTYPE(1:10)
      CHARACTER*4       ICRSH(0:2)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          A00AAF, DCOPY, E04NCN, F06DFF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NC/NEWOPT
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (MSGLS,MSGLVL)
C     .. Save statement ..
      SAVE              /AX02ZA/, /BE04NC/, /DE04NC/, /EE04NC/
C     .. Data statements ..
      DATA              CHESS(0), CHESS(1)/' NO', 'YES'/
      DATA              ICRSH(0), ICRSH(1), ICRSH(2)/'COLD', 'WARM',
     *                  'HOT '/
      DATA              LSTYPE(1), LSTYPE(2)/' FP', ' LP'/
      DATA              LSTYPE(3), LSTYPE(4), LSTYPE(5),
     *                  LSTYPE(6)/'QP1', 'QP2', 'QP3', 'QP4'/
      DATA              LSTYPE(7), LSTYPE(8), LSTYPE(9),
     *                  LSTYPE(10)/'LS1', 'LS2', 'LS3', 'LS4'/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
C
C     Make a dummy call to E04NCN to ensure that the defaults are set.
C
      CALL E04NCN(NOUT,'*',KEY)
      NEWOPT = .TRUE.
C
C     Save the optional parameters set by the user.  The values in
C     RPRMLS and IPRMLS may be changed to their default values.
C
      CALL F06DFF(MXPARM,IPRMLS,1,IPSVLS,1)
      CALL DCOPY(MXPARM,RPRMLS,1,RPSVLS,1)
C
      IF (MSGLVL.EQ.IDUMMY) MSGLVL = 10
      IF (IPRNT.LT.0) IPRNT = NOUT
      IF (ISUMRY.LT.0 .OR. MSGLVL.LT.5) ISUMRY = -1
      IPRINT = IPRNT
      ISUMM = ISUMRY
      IF (LPROB.LT.0) LPROB = 7
      IF (LCRASH.LT.0 .OR. LCRASH.GT.2) LCRASH = 0
      IF (LFORMH.LT.0 .OR. LFORMH.GT.1) LFORMH = 0
      IF (ITMAX1.LT.0) ITMAX1 = MAX(50,5*(N+NCLIN))
      IF (ITMAX2.LT.0) ITMAX2 = MAX(50,5*(N+NCLIN))
      IF (TOLACT.LT.ZERO) TOLACT = WRKTOL
      IF (TOLFEA.EQ.RDUMMY .OR. (TOLFEA.GE.ZERO .AND. TOLFEA.LT.EPSMCH))
     *    TOLFEA = EPSPT5
      IF (TOLRNK.LE.ZERO .AND. (LPROB.EQ.5 .OR. LPROB.EQ.7 .OR.
     *    LPROB.EQ.9)) TOLRNK = HUNDRD*EPSMCH
      IF (TOLRNK.LE.ZERO) TOLRNK = TEN*EPSPT5
      IF (BIGBND.LE.ZERO) BIGBND = GIGANT
      IF (BIGDX.LE.ZERO) BIGDX = MAX(GIGANT,BIGBND)
C
      IF (MSGLVL.GT.0) THEN
C
C        Print the title.
C
         LENT = LEN(TITLE)
         WRITE (REC,FMT='(/80A1)') (TITLE(J:J),J=1,LENT)
         CALL X04BAY(IPRINT,2,REC)
         CALL A00AAF
C
         IF (MSGLVL.GE.5 .AND. ISUMM.GE.0 .AND. ISUMM.NE.IPRINT) THEN
            WRITE (REC,FMT=99994) (TITLE(J:J),J=1,LENT)
            CALL X04BAY(ISUMM,2,REC)
         END IF
C
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99998) LSTYPE(LPROB), CHESS(LFORMH)
         CALL X04BAY(IPRINT,2,REC)
         WRITE (REC,FMT=99997) NCLIN, TOLFEA, N, TOLACT, M, TOLRNK
         CALL X04BAY(IPRINT,4,REC)
         WRITE (REC,FMT=99996) BIGBND, ICRSH(LCRASH), BIGDX, EPSMCH
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99995) MSGLVL, ITMAX1, ISUMM, ITMAX2
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      RETURN
C
C
C     End of  E04NCS. (LSDFLT)
C
99999 FORMAT (/' Parameters',/' ----------')
99998 FORMAT (/' Problem type...........',7X,A3,7X,'Hessian...........',
     *       '.....',7X,A3)
99997 FORMAT (/' Linear constraints.....',I10,7X,'Feasibility toleranc',
     *       'e..',1P,D10.2,6X,/' Variables..............',I10,7X,'Cra',
     *       'sh tolerance........',1P,D10.2,
     *       /' Objective matrix rows..',I10,7X,'Rank tolerance.......',
     *       '..',1P,D10.2)
99996 FORMAT (/' Infinite bound size....',1P,D10.2,7X,A4,' start......',
     *       '.......',/' Infinite step size.....',1P,D10.2,7X,'EPS (m',
     *       'achine precision)',1P,D10.2)
99995 FORMAT (/' Print level............',I10,7X,'Feasibility phase it',
     *       'ns.',I10,/' Monitoring file........',I10,7X,'Optimality ',
     *       ' phase itns.',I10)
99994 FORMAT (/11A1,' monitoring information ')
      END
