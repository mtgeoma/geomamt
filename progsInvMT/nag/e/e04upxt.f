      SUBROUTINE E04UPX(M,N,NCLIN,NCNLN,TITLE)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1103 (JUL 1993).
C     MARK 17 REVISED. IER-1623 (JUN 1995).
C
C     ******************************************************************
C     E04UPX  loads the default values of parameters not set in the
C     options file.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written 10-September-1985.
C     This version of E04UPX dated  12-Jul-94.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  POINT3, POINT8
      PARAMETER         (POINT3=3.3D-1,POINT8=0.8D+0)
      DOUBLE PRECISION  POINT9, TWO
      PARAMETER         (POINT9=0.9D+0,TWO=2.0D+0)
      DOUBLE PRECISION  TENP6, HUNDRD
      PARAMETER         (TENP6=1.0D+6,HUNDRD=10.0D+1)
      DOUBLE PRECISION  RDUMMY
      INTEGER           IDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0,IDUMMY=-11111)
      DOUBLE PRECISION  GIGANT
      PARAMETER         (GIGANT=1.0D+20*0.99999D+0)
      DOUBLE PRECISION  WRKTOL
      PARAMETER         (WRKTOL=1.0D-2)
C     .. Scalar Arguments ..
      INTEGER           M, N, NCLIN, NCNLN
      CHARACTER*(*)     TITLE
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT, CTOL,
     *                  DXLIM, EPSPT3, EPSPT5, EPSPT8, EPSPT9, EPSRF,
     *                  ETA, FDINT, FTOL, HCNDBD, TOLACT, TOLFEA, TOLRNK
      INTEGER           IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1, ITMAX2,
     *                  ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4, KSAVE,
     *                  LCRASH, LFDSET, LFORMH, LINES1, LINES2, LPROB,
     *                  LTYPEH, LVERFY, LVLDER, LVLDIF, LVRFYC, MSGLS,
     *                  MSGNP, NCDIFF, NFDIFF, NLNF, NLNJ, NLNX, NLOAD,
     *                  NN, NNCLIN, NNCNLN, NOUT, NPROB, NRESET, NSAVE
      LOGICAL           NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNL(30), RPADNP(22),
     *                  RPSVLS(MXPARM), RPSVNL(MXPARM), RPSVNP(MXPARM),
     *                  WMACH(15)
      INTEGER           IPADLS(19), IPADNL(28), IPADNP(15),
     *                  IPSVLS(MXPARM), IPSVNL(MXPARM), IPSVNP(MXPARM),
     *                  JVERFY(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  CONDBD, DCTOL, EPSMCH
      INTEGER           J, LENT, MSGQP, NCTOTL, NMAJOR, NMINOR, NPLIN
      CHARACTER*16      KEY
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNL(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNL(MXPARM), IPRMNP(MXPARM)
      CHARACTER*3       CHESS(0:1)
      CHARACTER*4       ICRSH(0:2)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          A00AAF, DCOPY, E04UPU, F06DFF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, LEN, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /BE04UP/IPSVNL, LTYPEH, NRESET, IPADNL
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /CE04UP/RPSVNL, RPADNL
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /EE04UC/NEWOPT
      COMMON            /GE04UC/IPSVNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, LVLDER, LVERFY, MSGNP, NLNF, NLNJ, NLNX,
     *                  NNCNLN, NSAVE, NLOAD, KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),ITMXNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (IPRMNL(1),LTYPEH), (RPRMNL(1),RPADNL(1))
      EQUIVALENCE       (ITMXNP,NMAJOR), (ITMAX2,NMINOR), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /AX02ZA/, /EE04UC/, /DE04NC/, /EE04NC/,
     *                  /GE04UC/, /HE04UC/, /BE04UP/, /CE04UP/
C     .. Data statements ..
      DATA              ICRSH(0), ICRSH(1), ICRSH(2)/'COLD', 'WARM',
     *                  'HOT '/
      DATA              CHESS(0), CHESS(1)/' NO', 'YES'/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
      NOUT = WMACH(11)
C
      CONDBD = MAX(ONE/(HUNDRD*EPSMCH*DBLE(N)),TENP6)
C
      NPLIN = N + NCLIN
      NCTOTL = NPLIN + NCNLN
C
C     Make a dummy call E04UPU to ensure that the defaults are set.
C
      CALL E04UPU(NOUT,'*',KEY)
      NEWOPT = .TRUE.
C
C     Save the optional parameters set by the user.  The values in
C     IPRMLS, RPRMLS, IPRMNP, RPRMNP, IPRMNL and RPRMNL may be changed
C     to their default values.
C
      CALL F06DFF(MXPARM,IPRMLS,1,IPSVLS,1)
      CALL DCOPY(MXPARM,RPRMLS,1,RPSVLS,1)
      CALL F06DFF(MXPARM,IPRMNP,1,IPSVNP,1)
      CALL DCOPY(MXPARM,RPRMNP,1,RPSVNP,1)
      CALL F06DFF(MXPARM,IPRMNL,1,IPSVNL,1)
      CALL DCOPY(MXPARM,RPRMNL,1,RPSVNL,1)
C
      IF (MSGNP.EQ.IDUMMY) MSGNP = 10
      IF (MSGQP.EQ.IDUMMY) MSGQP = 0
      IF (IPRNT.LT.0) IPRNT = NOUT
      IF (ISUMRY.LT.0 .OR. (MSGNP.LT.5 .AND. MSGQP.LT.5)) ISUMRY = -1
      IPRINT = IPRNT
      ISUMM = ISUMRY
      IF (LCRASH.LT.0 .OR. LCRASH.GT.2) LCRASH = 0
      IF (LVLDER.LT.0 .OR. LVLDER.GT.3) LVLDER = 3
      IF (LFORMH.LT.0 .OR. LFORMH.GT.1) LFORMH = 0
      IF (NMAJOR.LT.0) NMAJOR = MAX(50,3*NPLIN+10*NCNLN)
      IF (NMINOR.LT.1) NMINOR = MAX(50,3*NCTOTL)
      NLNF = N
      NLNJ = N
      NLNX = N
      IF (JVRFY2.LE.0 .OR. JVRFY2.GT.N) JVRFY2 = N
      IF (JVRFY1.LE.0 .OR. JVRFY1.GT.JVRFY2) JVRFY1 = 1
      IF (JVRFY4.LE.0 .OR. JVRFY4.GT.N) JVRFY4 = N
      IF (JVRFY3.LE.0 .OR. JVRFY3.GT.JVRFY4) JVRFY3 = 1
      IF ((LVERFY.LT.-1 .OR. LVERFY.GT.13)
     *    .OR. (LVERFY.GE.4 .AND. LVERFY.LE.9)) LVERFY = 0
C
      IF (TOLACT.LT.ZERO .OR. TOLACT.GE.ONE) TOLACT = WRKTOL
      IF (TOLFEA.LT.EPSMCH .OR. TOLFEA.GE.ONE) TOLFEA = EPSPT5
      IF (EPSRF.LT.EPSMCH .OR. EPSRF.GE.ONE) EPSRF = EPSPT9
      LFDSET = 0
      IF (FDINT.LT.ZERO) LFDSET = 2
      IF (FDINT.EQ.RDUMMY) LFDSET = 0
      IF (FDINT.GE.EPSMCH .AND. FDINT.LT.ONE) LFDSET = 1
      IF (LFDSET.EQ.1 .AND. (CDINT.LT.EPSMCH .OR. CDINT.GE.ONE))
     *    CDINT = EPSRF**POINT3
      IF (BIGBND.LE.ZERO) BIGBND = GIGANT
      IF (BIGDX.LE.ZERO) BIGDX = MAX(GIGANT,BIGBND)
      IF (DXLIM.LE.ZERO) DXLIM = TWO
      IF (ETA.LT.ZERO .OR. ETA.GE.ONE) ETA = POINT9
      IF (FTOL.LT.EPSRF .OR. FTOL.GE.ONE) FTOL = EPSRF**POINT8
C
      IF (HCNDBD.LT.ONE) HCNDBD = CONDBD
C
      DCTOL = EPSPT5
      IF (LVLDER.LT.2) DCTOL = EPSPT3
      IF (CTOL.LT.EPSMCH .OR. CTOL.GE.ONE) CTOL = DCTOL
      IF (LTYPEH.LT.0) LTYPEH = 0
      IF (NRESET.LE.0) NRESET = 2
C
      ITMAX1 = MAX(50,3*(N+NCLIN+NCNLN))
      JVERFY(1) = JVRFY1
      JVERFY(2) = JVRFY2
      JVERFY(3) = JVRFY3
      JVERFY(4) = JVRFY4
C
      IF (MSGNP.GT.0) THEN
C
C        Print the title.
C
         LENT = LEN(TITLE)
         WRITE (REC,FMT=99986) (TITLE(J:J),J=1,LENT)
         CALL X04BAY(IPRINT,2,REC)
         CALL A00AAF
C
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99998) NCLIN, N, NCNLN, M
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99997) BIGBND, ICRSH(LCRASH), BIGDX, EPSMCH,
     *     DXLIM, CHESS(LFORMH)
         CALL X04BAY(IPRINT,4,REC)
         WRITE (REC,FMT=99996) TOLFEA, TOLACT, CTOL, FTOL, ETA, EPSRF
         CALL X04BAY(IPRINT,4,REC)
         WRITE (REC,FMT=99995) LVLDER, ISUMRY, LVERFY
         CALL X04BAY(IPRINT,3,REC)
         IF (LVERFY.GT.0) THEN
            WRITE (REC,FMT=99994) JVRFY1, JVRFY2
            CALL X04BAY(IPRINT,2,REC)
            IF (NCNLN.GT.0) THEN
               WRITE (REC,FMT=99993) JVRFY3, JVRFY4
               CALL X04BAF(IPRINT,REC(1))
            END IF
         END IF
         WRITE (REC,FMT=99992) NMAJOR, MSGNP, NMINOR, MSGQP
         CALL X04BAY(IPRINT,3,REC)
C
         IF (LVLDER.LT.3) THEN
            IF (LFDSET.EQ.0) THEN
               WRITE (REC,FMT=99991)
               CALL X04BAY(IPRINT,2,REC)
            ELSE IF (LFDSET.EQ.1) THEN
               WRITE (REC,FMT=99990) FDINT, CDINT
               CALL X04BAY(IPRINT,2,REC)
            ELSE IF (LFDSET.EQ.2) THEN
               WRITE (REC,FMT=99989)
               CALL X04BAY(IPRINT,2,REC)
            END IF
         END IF
C
         IF (LTYPEH.EQ.0) THEN
            WRITE (REC,FMT=99988) NRESET
            CALL X04BAY(IPRINT,2,REC)
         ELSE
            WRITE (REC,FMT=99987) NRESET
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
      END IF
C
      IF (MSGNP.GE.5 .OR. MSGQP.GE.5) THEN
         IF (ISUMM.GE.0 .AND. ISUMM.NE.IPRINT) THEN
            LENT = LEN(TITLE)
            WRITE (REC,FMT=99985) (TITLE(J:J),J=1,LENT)
            CALL X04BAY(ISUMM,2,REC)
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04UPX.  (NLDFLT)
C
99999 FORMAT (/' Parameters',/' ----------')
99998 FORMAT (/' Linear constraints.....',I10,7X,'Variables...........',
     *       '...',I10,/' Nonlinear constraints..',I10,7X,'Subfunction',
     *       's...........',I10)
99997 FORMAT (/' Infinite bound size....',1P,D10.2,7X,A4,' start......',
     *       '.......',/' Infinite step size.....',1P,D10.2,7X,'EPS (m',
     *       'achine precision)',1P,D10.2,/' Step limit.............',
     *       1P,D10.2,7X,'Hessian................',7X,A3)
99996 FORMAT (/' Linear feasibility.....',1P,D10.2,7X,'Crash tolerance',
     *       '........',1P,D10.2,/' Nonlinear feasibility..',1P,D10.2,
     *       7X,'Optimality tolerance...',1P,D10.2,/' Line search tole',
     *       'rance..',1P,D10.2,7X,'Function precision.....',1P,D10.2)
99995 FORMAT (/' Derivative level.......',I10,7X,'Monitoring file.....',
     *       '...',I10,/' Verify level...........',I10)
99994 FORMAT (/' Start obj chk at varble',I10,7X,'Stop obj chk at varb',
     *       'le.',I10)
99993 FORMAT (' Start con chk at varble',I10,7X,'Stop con chk at varbl',
     *       'e.',I10)
99992 FORMAT (/' Major iterations limit.',I10,7X,'Major print level...',
     *       '...',I10,/' Minor iterations limit.',I10,7X,'Minor print',
     *       ' level......',I10)
99991 FORMAT (/' Difference intervals to be computed.')
99990 FORMAT (/' Difference interval....',1P,D10.2,7X,'Central diffce ',
     *       'interval',1P,D10.2)
99989 FORMAT (/' User-supplied difference intervals.')
99988 FORMAT (/' J''J initial Hessian....',17X,
     *       'Reset frequency........',I10)
99987 FORMAT (/' Unit initial Hessian...',17X,'Reset frequency........',
     *       I10)
99986 FORMAT (/80A1)
99985 FORMAT (/11A1,' monitoring information ')
      END
