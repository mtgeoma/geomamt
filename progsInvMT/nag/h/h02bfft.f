      SUBROUTINE H02BFF(INFILE,MAXN,MAXM,OPTIM,XBLDEF,XBUDEF,MAXDPT,
     *                  MSGLVL,N,M,X,CRNAME,IWORK,LIWORK,RWORK,LRWORK,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='H02BFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XBLDEF, XBUDEF
      INTEGER           IFAIL, INFILE, LIWORK, LRWORK, M, MAXDPT, MAXM,
     *                  MAXN, MSGLVL, N
      CHARACTER*3       OPTIM
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), X(MAXN)
      INTEGER           IWORK(LIWORK)
      CHARACTER*8       CRNAME(MAXN+MAXM)
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, TOLACT, TOLFEA,
     *                  TOLRNK
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(30)
C     .. Local Scalars ..
      DOUBLE PRECISION  INFNTY, OBJMIP, OBJVAL, TOLFES, TOLIV
      INTEGER           I, IAPTR, IAX, IBLPTR, IBUPTR, ICLPTR, ICVPTR,
     *                  IFAILX, IIDIM, IINPTR, IIPBU, IIPCV, IIPTR,
     *                  IISPTR, INTFST, IPBL, IPBU, IPCV, IPIS, IRDIM,
     *                  IREC, IRPTR, ISIZE, ITER, ITMAX, LDA, MAXNOD,
     *                  MSG, MSGLV1, NCTOTL, NOUT
      LOGICAL           MPSLST
      CHARACTER*3       LOPTIM
      CHARACTER*8       KBLANK, NMBND, NMOBJ, NMPROB, NMRHS, NMRNG
      CHARACTER*16      STR
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04MFF, E04MHF, E04UDU, H02BBF, H02BFZ, H02BUF,
     *                  H02BVF, X04ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /GE04MF/RPSVLC, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLC
C     .. Save statement ..
      SAVE              /GE04MF/
C     .. Data statements ..
      DATA              KBLANK/'        '/
C     .. Executable Statements ..
C
C     INITIALIZE PARAMETERS
      LOPTIM = OPTIM
C     Convert to upper case
      CALL E04UDU(LOPTIM)
C
      N = MAXN
      M = MAXM
      IREC = 0
      IFAILX = 0
      IF (INFILE.LT.0 .OR. INFILE.GT.99) THEN
         IFAILX = 5
         WRITE (REC,FMT=99999) INFILE
         IREC = 3
         GO TO 80
      END IF
      IF (MAXN.LT.1) THEN
         IFAILX = 5
         WRITE (REC,FMT=99998) MAXN
         IREC = 3
         GO TO 80
      END IF
      IF (MAXM.LT.1) THEN
         IFAILX = 5
         WRITE (REC,FMT=99997) MAXM
         IREC = 3
         GO TO 80
      END IF
      IF (XBLDEF.GT.XBUDEF) THEN
         IFAILX = 5
         WRITE (REC,FMT=99996) XBLDEF, XBUDEF
         IREC = 3
         GO TO 80
      END IF
      IF (LOPTIM.NE.'MIN' .AND. LOPTIM.NE.'MAX') THEN
         IFAILX = 5
         WRITE (REC,FMT=99995) OPTIM
         IREC = 3
         GO TO 80
      END IF
      CALL X04ABF(0,NOUT)
      NMPROB = KBLANK
      NMOBJ = KBLANK
      NMRHS = KBLANK
      NMRNG = KBLANK
      NMBND = KBLANK
      LDA = M
      MPSLST = .FALSE.
C
      CALL E04MHF('Nolist')
      CALL E04MHF('Infinite Bound Size = 1.0D+20   * (i.e. BIGBND) ')
C
C     SET UP POINTERS (REAL)
C
      IAPTR = 1
      IBLPTR = IAPTR + M*N
      IBUPTR = IBLPTR + N + M
      ICLPTR = IBUPTR + N + M
      ICVPTR = ICLPTR + N + M
      IAX = ICVPTR + N
      IRPTR = IAX + N + M
C
C     SET UP POINTERS (INTEGER)
C
      IINPTR = 1
      IISPTR = IINPTR + N
      IIPTR = IISPTR + N + M
C
      IIDIM = LIWORK - IIPTR + 1
      IF (IIDIM.LT.1) THEN
         IFAILX = 5
         WRITE (REC,FMT=99994) LIWORK
         IREC = 3
         GO TO 80
      END IF
C
      IRDIM = LRWORK - IRPTR + 1
      IF (IRDIM.LT.1) THEN
         IFAILX = 5
         WRITE (REC,FMT=99993) LRWORK
         IREC = 3
         GO TO 80
      END IF
C
      IFAILX = -13
      IF (IFAIL.EQ.1) IFAILX = 1
C
      CALL H02BUF(INFILE,MAXN,MAXM,OPTIM,XBLDEF,XBUDEF,NMOBJ,NMRHS,
     *            NMRNG,NMBND,MPSLST,N,M,RWORK(IAPTR),RWORK(IBLPTR),
     *            RWORK(IBUPTR),RWORK(ICVPTR),X,IWORK(IINPTR),CRNAME,
     *            NMPROB,IWORK(IISPTR),IFAILX)
C
      IFAILX = -IFAILX
      IF (IFAILX.NE.0) GO TO 80
C
C     DETERMINE WHICH SOLVER TO USE
C
      DO 20 I = 1, N
         IF (IWORK(IINPTR+I-1).EQ.1) GO TO 40
   20 CONTINUE
C
C     SOLVE THE LP PROBLEM USING E04MFF
C
      ISIZE = 4*MAXN + MAXM + 3
      IF (LIWORK.LT.ISIZE) THEN
         IFAILX = 5
         IREC = 4
         WRITE (REC,FMT=99991) LIWORK, ISIZE
         GO TO 80
      END IF
      ISIZE = 2*MIN(MAXN,MAXM+1)**2 + MAXM*MAXN + 12*MAXN + 9*MAXM
      IF (LRWORK.LT.ISIZE) THEN
         IFAILX = 5
         IREC = 4
         WRITE (REC,FMT=99990) LRWORK, ISIZE
         GO TO 80
      END IF
C
      WRITE (STR,FMT='(I5)') MSGLVL
      CALL E04MHF('Print Level = '//STR)
C     CHANGE MSGLVL
      IF (MSGLVL.LT.2) THEN
         MSG = 0
         WRITE (STR,FMT='(I5)') MSG
         CALL E04MHF('Print Level = '//STR)
      ELSE IF (MSGLVL.EQ.10) THEN
         MSG = 5
         WRITE (STR,FMT='(I5)') MSG
         CALL E04MHF('Print Level = '//STR)
      END IF
      IFAILX = -13
      IF (IFAIL.EQ.1) IFAILX = 1
C
      CALL E04MFF(N,M,RWORK(IAPTR),LDA,RWORK(IBLPTR),RWORK(IBUPTR),
     *            RWORK(ICVPTR),IWORK(IISPTR),X,ITER,OBJVAL,RWORK(IAX),
     *            RWORK(ICLPTR),IWORK(IIPTR),IIDIM,RWORK(IRPTR),IRDIM,
     *            IFAILX)
C
      IF (IFAILX.EQ.6) GO TO 80
C
      IF (MSGLVL.GT.0 .AND. MSGLVL.NE.5) CALL H02BVF(N,M,RWORK(IAPTR),
     *    LDA,RWORK(IBLPTR),RWORK(IBUPTR),X,RWORK(ICLPTR),IWORK(IISPTR),
     *    CRNAME,IFAIL)
C
C     SAVE THE SOLUTION
C
      NCTOTL = N + M
      CALL H02BFZ(RWORK(IBLPTR),RWORK(IBUPTR),RWORK(ICLPTR),
     *            IWORK(IISPTR),NCTOTL,IWORK,LIWORK,RWORK,LRWORK)
      GO TO 80
C
   40 IF (MAXDPT.LT.2) THEN
         IFAILX = 5
         IREC = 3
         WRITE (REC,FMT=99992) MAXDPT
         GO TO 80
      END IF
      ISIZE = (25+MAXN+MAXM)*MAXDPT + 2*MAXM + 7*MAXN + 4
      IF (LIWORK.LT.ISIZE) THEN
         IFAILX = 5
         IREC = 4
         WRITE (REC,FMT=99991) LIWORK, ISIZE
         GO TO 80
      END IF
      ISIZE = MAXDPT*(MAXN+2) + 2*MIN(MAXN,MAXM+1)**2 + MAXM*MAXN +
     *        18*MAXN + 15*MAXM
      IF (LRWORK.LT.ISIZE) THEN
         IFAILX = 5
         IREC = 4
         WRITE (REC,FMT=99990) LRWORK, ISIZE
         GO TO 80
      END IF
      ITMAX = -11111
      MAXNOD = 0
      INTFST = 0
      TOLFES = 0.0D+0
      TOLIV = 1.0D-5
      MSGLV1 = MSGLVL
C
C     SOLVE THE IP PROBLEM USING H02BBF
C
      IFAILX = -13
      IF (IFAIL.EQ.1) IFAILX = 1
C
      INFNTY = BIGBND
      CALL H02BBF(ITMAX,MSGLV1,N,M,RWORK(IAPTR),LDA,RWORK(IBLPTR),
     *            RWORK(IBUPTR),IWORK(IINPTR),RWORK(ICVPTR),MAXNOD,
     *            INTFST,MAXDPT,TOLIV,TOLFES,INFNTY,X,OBJMIP,
     *            IWORK(IIPTR),IIDIM,RWORK(IRPTR),IRDIM,IFAILX)
C
      IF (IFAILX.NE.0) THEN
         IF (IFAILX.EQ.1) IFAILX = 9
         IF (IFAILX.NE.7) GO TO 80
      END IF
C
      NCTOTL = N + M
C
      IPBL = IRPTR
      IPBU = IRPTR + NCTOTL
      IPCV = IPBU + NCTOTL
      IPIS = IIPTR
      IF (MSGLVL.GT.0) CALL H02BVF(N,M,RWORK(IAPTR),LDA,RWORK(IPBL),
     *                             RWORK(IPBU),X,RWORK(IPCV),IWORK(IPIS)
     *                             ,CRNAME,IFAIL)
C
C     SAVE THE INTEGER SOLUTION
C
      IIPBU = NCTOTL
      IIPCV = IIPBU + NCTOTL
      DO 60 I = 1, NCTOTL
         IWORK(I) = IWORK(IPIS+I-1)
         RWORK(I) = RWORK(IPBL+I-1)
         RWORK(I+IIPBU) = RWORK(IPBU+I-1)
         RWORK(I+IIPCV) = RWORK(IPCV+I-1)
   60 CONTINUE
C
   80 IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,IREC,REC)
C
      RETURN
C
99999 FORMAT (/' ** On entry, INFILE is out of range:',/'    INFILE = ',
     *       I16)
99998 FORMAT (/' ** On entry, MAXN.lt.1:',/'    MAXN = ',I16)
99997 FORMAT (/' ** On entry, MAXM.lt.1:',/'    MAXM = ',I16)
99996 FORMAT (/' ** On entry, XBLDEF.gt.XBUDEF:',/'    XBLDEF = ',G16.7,
     *       '   XBUDEF = ',G16.7)
99995 FORMAT (/' ** On entry, OPTIM does not equal MIN or MAX:',/'    ',
     *       'OPTIM = ',A3)
99994 FORMAT (/' ** On entry, not enough integer workspace to read dat',
     *       'a file:',/'    LIWORK = ',I16)
99993 FORMAT (/' ** On entry, not enough real workspace to read data f',
     *       'ile:',/'    LRWORK = ',I16)
99992 FORMAT (/' ** On entry, MAXDPT.lt.2:',/'    MAXDPT = ',I16)
99991 FORMAT (/' ** On entry, not enough integer workspace to solve pr',
     *       'oblem:',/'    LIWORK = ',I16,/'    LIWORK must be at lea',
     *       'st ',I16)
99990 FORMAT (/' ** On entry, not enough real workspace to solve probl',
     *       'em:',/'    LRWORK = ',I16,/'    LRWORK must be at least ',
     *       I16)
      END
