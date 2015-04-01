      SUBROUTINE F11GAF(METHOD,PRECON,SIGCMP,NORM,WEIGHT,ITERM,N,TOL,
     *                  MAXITN,ANORM,SIGMAX,SIGTOL,MAXITS,MONIT,LWREQ,
     *                  IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GAF - Initialization for the symmetric iterative solver suite
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F11GAF')
      DOUBLE PRECISION  ZERO, ONE, HUNDTH
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HUNDTH=1.0D-2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANORM, SIGMAX, SIGTOL, TOL
      INTEGER           IFAIL, ITERM, LWREQ, MAXITN, MAXITS, MONIT, N
      CHARACTER         NORM, PRECON, SIGCMP, WEIGHT
      CHARACTER*(*)     METHOD
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORL, SIGML, SIGTLL, TOLL
      INTEGER           IMETH, INFO, INORM, IPREC, ISIGC, IWEIG, MAXITL,
     *                  NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  RDATA(20)
      INTEGER           IDATA(20)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F06FBF, F11BAZ
C     .. Intrinsic Functions ..
      INTRINSIC         INDEX, MAX, SIGN, SQRT
C     .. Executable Statements ..
C
C     Initialize
C
      INFO = 0
      NREC = 0
      LWREQ = 0
      IF (IFAIL.NE.0) IFAIL = SIGN(1,IFAIL)
C
C     Check the arguments
C
      IF (INDEX(METHOD,'CG').EQ.1 .OR. INDEX(METHOD,'cg').EQ.1) THEN
         IMETH = 1
      ELSE IF (INDEX(METHOD,'SYMMLQ').EQ.1 .OR. INDEX(METHOD,'symmlq')
     *         .EQ.1) THEN
         IMETH = 2
      ELSE
         IMETH = -1
      END IF
      IPREC = (INDEX('NnPp',PRECON)+1)/2 - 1
      ISIGC = (INDEX('NnSs',SIGCMP)+1)/2 - 1
      INORM = (INDEX('11Ii22',NORM)+1)/2 - 1
      IWEIG = (INDEX('NnWw',WEIGHT)+1)/2 - 1
      TOLL = MAX(TOL,ZERO)
C
      IF (IMETH.LT.0) THEN
         INFO = -1
      ELSE IF (IPREC.LT.0) THEN
         INFO = -2
      ELSE IF (ISIGC.LT.0) THEN
         INFO = -3
      ELSE IF (INORM.LT.0) THEN
         INFO = -4
      ELSE IF (IWEIG.LT.0) THEN
         INFO = -5
C      ELSE IF ((ITERM.NE.1) .AND. (ITERM.NE.2)) THEN
      ELSE IF ((ITERM.NE.1) .AND. ((ITERM.NE.2) .OR. (IMETH.LE.1))) THEN
         INFO = -6
      ELSE IF (N.LE.0) THEN
         INFO = -7
      ELSE IF (TOLL.GE.ONE) THEN
         INFO = -8
      ELSE IF (MAXITN.LE.0) THEN
         INFO = -9
      ELSE IF (MONIT.GT.MAXITN) THEN
         INFO = -14
      END IF
      IF (INFO.NE.0) GO TO 20
C
      IF (ISIGC.GE.1) THEN
         ISIGC = 2
      ELSE IF ((ITERM.GE.2) .OR. ((IMETH.GE.2) .AND. (IWEIG.LE.0))) THEN
         ISIGC = 1
      END IF
      IF (ITERM.NE.2) THEN
         ANORL = MAX(ANORM,ZERO)
      ELSE
         ANORL = ZERO
      END IF
      IF ((ITERM .NE. 2) .AND. (ISIGC.LE.0)) THEN
         SIGML = ZERO
      ELSE
         SIGML = MAX(SIGMAX,ZERO)
      END IF
      IF (ISIGC.LE.0) THEN
         SIGTLL = ZERO
         MAXITL = 0
      ELSE
         IF (SIGML.GT.ZERO) ISIGC = 0
         IF (ISIGC.LE.1) THEN
            SIGTLL = ZERO
            MAXITL = 0
         ELSE
            SIGTLL = MAX(SIGTOL,ZERO)
            MAXITL = MAXITS
         END IF
      END IF
C
      IF ((ITERM.GE.2) .AND. (IWEIG.GE.1)) THEN
         INFO = -5
      ELSE IF ((ITERM.GE.2) .AND. (INORM.LE.1)) THEN
         INFO = -6
      ELSE IF ((ITERM.LE.1) .AND. (INORM.GE.2) .AND. (ANORL.LE.ZERO))
     *         THEN
         INFO = -10
      ELSE IF ((ISIGC.GE.2) .AND. (SIGMAX.LE.ZERO)) THEN
         IF (SIGTLL.GE.ONE) THEN
            INFO = -12
         ELSE IF ((MAXITL.LE.0) .OR. (MAXITL.GT.MAXITN)) THEN
            INFO = -13
         END IF
      END IF
      IF (INFO.NE.0) GO TO 20
C
C     Complete and store the parameters
C
      IF (TOLL.LE.ZERO) THEN
         TOLL = SQRT(X02AJF())
      ELSE
         TOLL = MAX(X02AJF(),TOLL)
      END IF
C
      IF (ISIGC.GE.2) THEN
         IF (SIGTLL.LE.ZERO) THEN
            SIGTLL = HUNDTH
         ELSE
            SIGTLL = MAX(X02AJF(),SIGTLL)
         END IF
      END IF
C
      IF (IMETH.LE.1) THEN
         LWREQ = 5*N
      ELSE
         LWREQ = 6*N
      END IF
      IF (IWEIG.GE.1) LWREQ = LWREQ + N
      IF (ISIGC.GE.2) LWREQ = LWREQ + 2*(MAXITL+1)
C
      CALL F06DBF(20,0,IDATA,1)
      CALL F06FBF(20,ZERO,RDATA,1)
      IDATA(2) = IMETH
      IDATA(3) = IPREC
      IDATA(4) = ISIGC
      IDATA(5) = INORM + 3*IWEIG
      IDATA(6) = ITERM
      IDATA(7) = N
      IDATA(8) = MAXITN
      IDATA(9) = MAXITL
      IDATA(10) = MAX(MONIT,0)
      IDATA(11) = LWREQ
C
      RDATA(1) = TOLL
      RDATA(2) = SIGTLL
      RDATA(3) = ANORL
      RDATA(4) = SIGML
C
C     Store the parameters
C
      CALL F11BAZ(1,IDATA,RDATA,INFO)
C
C     Complete
C
   20 CONTINUE
      IF (INFO.LT.0) THEN
         NREC = 1
         WRITE (REC(1),FMT=99999) - INFO
      ELSE IF (INFO.GT.0) THEN
         NREC = 2
         WRITE (REC,FMT=99998)
      END IF
C
      IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,REC)
C
C     End of subroutine F11GAF
C
      RETURN
C
99999 FORMAT (' ** On entry, parameter number ',I2,' had an illegal va',
     *       'lue.')
99998 FORMAT (' ** F11GAF has been called out of sequence: either ',
     *       'F11GAF has been called',/' ** twice or F11GBF has not',
     *       ' terminated its current task.')
      END
