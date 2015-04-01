      SUBROUTINE E04DGR(N,TITLE)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-1054 (JUL 1993).
C
C     E04DGR  loads the default values of parameters not set by the
C     user.
C
C     -- Written on 6-June-1986.
C     Sven Hammarling and Janet Welding, NAG Central Office.
C
C     AE04DG and BE04DG  are the integer and real common blocks for the
C     optional parameters.  CE04DG is the common block to pass NEWOPT
C     between here and the user-callable option setting routines.
C
C     .. Parameters ..
      INTEGER           MXPARM, NIPARM, NRPARM
      PARAMETER         (MXPARM=30,NIPARM=8,NRPARM=5)
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      DOUBLE PRECISION  POINT8, POINT9
      PARAMETER         (POINT8=0.8D+0,POINT9=0.9D+0)
      DOUBLE PRECISION  RDUMMY
      INTEGER           IDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0,IDUMMY=-11111)
      DOUBLE PRECISION  GIGANT
      PARAMETER         (GIGANT=1.0D+20*0.99999D0)
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER*(*)     TITLE
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGDX, EPSPT3, EPSPT5, EPSPT8, EPSPT9, EPSRF,
     *                  ETA, FGUESS, FTOL
      INTEGER           IDBGCG, IPRINT, ISUMM, ITMAX, JVRFY1, JVRFY2,
     *                  LDBGCG, LINES1, LINES2, LVERFY, MSGCG, NN, NOUT
      LOGICAL           NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADCG(MXPARM-NRPARM), RPSVCG(MXPARM), WMACH(15)
      INTEGER           IPADCG(MXPARM-NIPARM), IPSVCG(MXPARM)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSMCH
      INTEGER           IDBG, J, LENT, MSGDBG, MSGLVL
      CHARACTER*16      KEY
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMCG(MXPARM)
      INTEGER           IPRMCG(MXPARM)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          A00AAF, DCOPY, E04DGS, F06DFF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MAX
C     .. Common blocks ..
      COMMON            /AE04DG/IPSVCG, IDBGCG, ITMAX, JVRFY1, JVRFY2,
     *                  LDBGCG, LVERFY, MSGCG, NN, IPADCG
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04DG/RPSVCG, BIGDX, EPSRF, ETA, FGUESS, FTOL,
     *                  RPADCG
      COMMON            /CE04DG/NEWOPT
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Equivalences ..
      EQUIVALENCE       (IPRMCG(1),IDBGCG), (RPRMCG(1),BIGDX)
      EQUIVALENCE       (MSGCG,MSGLVL), (IDBGCG,IDBG), (LDBGCG,MSGDBG)
C     .. Save statement ..
      SAVE              /AE04DG/, /BE04DG/, /CE04DG/, /AX02ZA/
C     .. Executable Statements ..
C
C     Make a dummy call to  E04DGS  to ensure that the defaults are set
C     and set  NEWOPT  as true so that the next call to an option
C     routine recognizes that  E04DGF  has been called.
C
      CALL E04DGS(NOUT,'*',KEY)
      NEWOPT = .TRUE.
C
C     Save the optional parameters set by the user.  The values in
C     RPRMCG and IPRMCG may be changed to their default values.
C
      CALL F06DFF(MXPARM,IPRMCG,1,IPSVCG,1)
      CALL DCOPY(MXPARM,RPRMCG,1,RPSVCG,1)
C
C     Set the defaults.
C
      EPSMCH = WMACH(3)
      IPRINT = NOUT
      IF (BIGDX.LE.ZERO) BIGDX = GIGANT
      IF (ITMAX.LT.0) ITMAX = MAX(50,5*N)
      IF (MSGLVL.EQ.IDUMMY) MSGLVL = 10
      IF (MSGDBG.LT.0) MSGDBG = 0
      IF (MSGDBG.EQ.0) IDBG = ITMAX + 1
      IF (IDBG.LT.0) IDBG = 0
      IF ((JVRFY2.LE.0) .OR. (JVRFY2.GT.N)) JVRFY2 = N
      IF ((JVRFY1.LE.0) .OR. (JVRFY1.GT.JVRFY2)) JVRFY1 = 1
      IF ((LVERFY.LT.(-1)) .OR. (LVERFY.GT.1)) LVERFY = 0
      IF ((EPSRF.LT.EPSMCH) .OR. (EPSRF.GE.ONE)) EPSRF = EPSPT9
      IF ((FTOL.LT.EPSRF) .OR. (FTOL.GE.ONE)) FTOL = EPSRF**POINT8
      IF ((ETA.LT.ZERO) .OR. (ETA.GE.ONE)) ETA = POINT9
C
      IF (MSGLVL.GT.0) THEN
C
C        Print the title.
C
         LENT = LEN(TITLE)
         IF (LENT.GT.0) THEN
            WRITE (REC,FMT='( / ( 80A1 ) )') (TITLE(J:J),J=1,LENT)
            CALL X04BAY(NOUT,2,REC)
         END IF
         CALL A00AAF
C
C        Print the optional parameter values.
C
         WRITE (REC,FMT=99999)
         CALL X04BAY(NOUT,3,REC)
         WRITE (REC,FMT=99998) N
         CALL X04BAY(NOUT,2,REC)
         WRITE (REC,FMT=99997) BIGDX, EPSMCH, FTOL, ETA
         CALL X04BAY(NOUT,3,REC)
         IF (FGUESS.EQ.RDUMMY) THEN
            WRITE (REC,FMT=99996) EPSRF, LVERFY
         ELSE
            WRITE (REC,FMT=99995) FGUESS, EPSRF, LVERFY
         END IF
         CALL X04BAY(NOUT,3,REC)
         IF (LVERFY.GT.0) THEN
            WRITE (REC,FMT=99993) JVRFY1, JVRFY2
            CALL X04BAY(NOUT,2,REC)
         END IF
         WRITE (REC,FMT=99994) ITMAX, MSGLVL
         CALL X04BAY(NOUT,2,REC)
         IF (MSGDBG.GT.0) THEN
            WRITE (REC,FMT=99992) MSGDBG, IDBG
            CALL X04BAF(NOUT,REC(1))
         END IF
      END IF
C
      RETURN
C
C     End of  E04DGR. (CGDFLT)
C
99999 FORMAT (/' Parameters',/' ----------')
99998 FORMAT (/' Variables..............',I10)
99997 FORMAT (/' Maximum step length....',1P,D10.2,7X,'EPS (machine pr',
     *       'ecision)',1P,D10.2,/' Optimality tolerance...',1P,D10.2,
     *       7X,'Linesearch tolerance...',1P,D10.2)
99996 FORMAT (/' Est. opt. function val.      None',7X,'Function preci',
     *       'sion.....',1P,D10.2,/' Verify level...........',I10)
99995 FORMAT (/' Est. opt. function val.',1P,D10.2,7X,'Function precis',
     *       'ion.....',1P,D10.2,/' Verify level...........',I10)
99994 FORMAT (/' Iteration limit........',I10,7X,'Print level.........',
     *       '...',I10)
99993 FORMAT (/' Start obj chk at varble',I10,7X,'Stop obj chk at varb',
     *       'le.',I10)
99992 FORMAT (' Debug level............',I10,7X,'Debug start at itn...',
     *       '..',I10)
      END
