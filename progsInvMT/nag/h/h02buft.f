      SUBROUTINE H02BUF(INFILE,MAXN,MAXM,OPTIM,XBLDEF,XBUDEF,NMOBJ,
     *                  NMRHS,NMRNG,NMBND,MPSLST,N,M,A,BL,BU,CVEC,X,
     *                  INTVAR,CRNAME,NMPROB,IWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='H02BUF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XBLDEF, XBUDEF
      INTEGER           IFAIL, INFILE, M, MAXM, MAXN, N
      LOGICAL           MPSLST
      CHARACTER*3       OPTIM
      CHARACTER*8       NMBND, NMOBJ, NMPROB, NMRHS, NMRNG
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MAXM,MAXN), BL(MAXN+MAXM), BU(MAXN+MAXM),
     *                  CVEC(MAXN), X(MAXN)
      INTEGER           INTVAR(MAXN), IWORK(MAXN+MAXM)
      CHARACTER*8       CRNAME(MAXN+MAXM)
C     .. Local Scalars ..
      DOUBLE PRECISION  ROPTDR
      INTEGER           I, IERR, IOBJ, IREC, J, NCTOTL, NOUT
      CHARACTER*3       LOPTIM
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04UDU, H02BUY, H02BUZ, X04ABF
C     .. Executable Statements ..
C
C     INITIALIZE PARAMETERS
      LOPTIM = OPTIM
C     Convert to upper case
      CALL E04UDU(LOPTIM)
C
C     FORM THE DATA FOR THE PROBLEM
C
      IF (INFILE.LT.0 .OR. INFILE.GT.99) THEN
         IERR = 14
         WRITE (REC,FMT=99999) INFILE
         IREC = 3
         GO TO 40
      END IF
      IF (MAXN.LT.1) THEN
         IERR = 14
         WRITE (REC,FMT=99998) MAXN
         IREC = 3
         GO TO 40
      END IF
      IF (MAXM.LT.1) THEN
         IERR = 14
         WRITE (REC,FMT=99997) MAXM
         IREC = 3
         GO TO 40
      END IF
      IF (XBLDEF.GT.XBUDEF) THEN
         IERR = 14
         WRITE (REC,FMT=99995) XBLDEF, XBUDEF
         IREC = 3
         GO TO 40
      END IF
      IF (LOPTIM.NE.'MAX' .AND. LOPTIM.NE.'MIN') THEN
         IERR = 14
         WRITE (REC,FMT=99996) OPTIM
         IREC = 3
         GO TO 40
      END IF
      IREC = 2
      CALL X04ABF(0,NOUT)
      DO 20 J = 1, MAXN
         INTVAR(J) = 0
   20 CONTINUE
C
      IERR = 0
      CALL H02BUY(INFILE,NOUT,MPSLST,MAXM,MAXN,MAXN+MAXM,N,M,NCTOTL,
     *            NMPROB,NMOBJ,NMRHS,NMRNG,NMBND,IOBJ,XBLDEF,XBUDEF,
     *            CRNAME,A,BL,BU,X,INTVAR,IWORK,REC,IERR)
C
   40 IFAIL = P01ABF(IFAIL,IERR,SRNAME,IREC,REC)
      IF (IFAIL.NE.0) RETURN
C
      IF (LOPTIM.EQ.'MIN') THEN
         ROPTDR = 1.0D0
      ELSE
         ROPTDR = -1.0D0
      END IF
C
      DO 60 I = 1, N
         CVEC(I) = ROPTDR*A(IOBJ,I)
   60 CONTINUE
C
C     TAKE THE OBJECTIVE ROW INFO OUT
C
      CALL H02BUZ(MAXM,A,M,N,IOBJ,BL,BU,CRNAME)
      M = M - 1
      RETURN
C
99999 FORMAT (/' ** On entry, INFILE is out of range:',/'    INFILE = ',
     *       I16)
99998 FORMAT (/' ** On entry, MAXN.lt.1:',/'    MAXN = ',I16)
99997 FORMAT (/' ** On entry, MAXM.lt.1:',/'    MAXM = ',I16)
99996 FORMAT (/' ** On entry, OPTIM does not equal MIN or MAX:',/'    ',
     *       'OPTIM = ',A3)
99995 FORMAT (/' ** On entry, XBLDEF.gt.XBUDEF:',/'    XBLDEF = ',G16.7,
     *       '   XBUDEF = ',G16.7)
      END
