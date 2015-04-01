      SUBROUTINE H02BBF(ITMAX,MSGLVL,N,M,A,LDA,BL,BU,INTVAR,CVEC,MAXNOD,
     *                  INTFST,MAXDPT,TOLIV,TOLFES,BIGBND,X,OBJMIP,
     *                  IWORK,LIWORK,RWORK,LRWORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 14B REVISED. IER-843 (MAR 1990).
C     MARK 15 REVISED. IER-927 (APR 1991).
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1667 (JUN 1995).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='H02BBF')
      DOUBLE PRECISION  RZERO
      PARAMETER         (RZERO=0.0D+0)
      INTEGER           IZERO
      PARAMETER         (IZERO=0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, OBJMIP, TOLFES, TOLIV
      INTEGER           IFAIL, INTFST, ITMAX, LDA, LIWORK, LRWORK, M,
     *                  MAXDPT, MAXNOD, MSGLVL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(M+N), BU(M+N), CVEC(N),
     *                  RWORK(LRWORK), X(N)
      INTEGER           INTVAR(N), IWORK(LIWORK)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2
      INTEGER           ACTNOD, AOPTVL, AX, BL1, BLO, BU1, BUP, BVINDX,
     *                  CLAM, CLAMDA, FATHER, FEASPT, I, IERR, IOPTCL,
     *                  ISTATE, ITMAX1, J, K, KIWORK, KLWORK, KRWORK, L,
     *                  LFTSON, LIWRK, LRWRK, MSG, MSGLV1, NCLIN,
     *                  NCTOTL, NERR, NERROR, NODEL, NODKNT, NON1, NON2,
     *                  NONOD, NXANOD, PARNOD, RITSON, SOLN, SOLVED,
     *                  SPLTVL, STACK
      LOGICAL           OK
      CHARACTER*10      STR
C     .. Local Arrays ..
      CHARACTER*5       ID(2)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          A00AAF, E04MHF, F06DBF, F06FBF, H02BBP, H02BBT,
     *                  H02BBZ, X02ZAZ, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ID(1), ID(2)/'Varbl', 'L Con'/
C     .. Executable Statements ..
C
C     RESET MSGLVL IF H02BBF IS CALLED FROM H02BFF
C     (AS INDICATED BY ITMAX = -11111).
C
      IF (ITMAX.EQ.-11111) THEN
         MSG = MSGLVL
         IF (MSGLVL.LT.5) THEN
            MSG = 0
         ELSE IF (MSGLVL.GT.5) THEN
            MSG = 20
         ELSE
            MSG = 15
         END IF
         MSGLV1 = MSGLVL
         ITMAX1 = ITMAX
      ELSE
         ITMAX1 = 0
      END IF
C
C     INITIALIZE VARIABLES.
C
      CALL X02ZAZ
      NOUT = WMACH(11)
      NERR = WMACH(12)
      IF (N.LE.0) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99992) N
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
      END IF
C
      IF (M.LT.0) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99991) M
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
      END IF
C
      IF (LDA.LT.MAX(1,M)) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99999) LDA, M
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
      END IF
C
      DO 20 I = 1, N
         IF (INTVAR(I).EQ.0 .OR. INTVAR(I).EQ.1) GO TO 20
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99998) I, I, INTVAR(I)
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
   20 CONTINUE
C
      NCLIN = M
      NCTOTL = N + NCLIN
      IF (MAXDPT.EQ.1 .AND. LRWORK.EQ.1 .AND. LIWORK.EQ.1) THEN
         KLWORK = 3*N/2
      ELSE
         KLWORK = MAXDPT
      END IF
C
      IF (MAXDPT.LT.2) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99990) MAXDPT
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
      END IF
C
C     COMPUTE THE AMOUNT OF REAL AND INTEGER WORKSPACE ACTUALLY REQUIRED
C
      KRWORK = KLWORK*(N+2) + 2*MIN(N,NCLIN+1)**2 + 13*N + 12*NCLIN
      KIWORK = (25+NCTOTL)*KLWORK + NCTOTL + 4*N + 4
C
      IF (KRWORK.GT.LRWORK) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99994) LRWORK, KRWORK
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
      END IF
C
      IF (KIWORK.GT.LIWORK) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99993) LIWORK, KIWORK
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
      END IF
C
      CALL E04MHF('Nolist')
      IF (ITMAX.LE.0) ITMAX = MAX(5*(NCLIN+N),50)
      WRITE (STR,FMT='(I10)') ITMAX
      CALL E04MHF('Iteration Limit = '//STR)
      IF (TOLIV.LE.RZERO) TOLIV = 1.0D-5
      IF (BIGBND.LE.RZERO) BIGBND = 1.0D+20
      IF (TOLFES.LE.RZERO) TOLFES = SQRT(X02AJF())
C
C     CHECK THE BOUNDS ON ALL VARIABLES AND CONSTRAINTS.
C
      NERROR = 0
      DO 40 J = 1, N + NCLIN
         B1 = BL(J)
         B2 = BU(J)
         OK = B1 .LT. B2 .OR. (B1.EQ.B2 .AND. ABS(B1).LT.BIGBND)
         IF ( .NOT. OK) THEN
            NERROR = NERROR + 1
            IF (J.GT.N) THEN
               K = J - N
               L = 2
            ELSE
               K = J
               L = 1
            END IF
            IF (B1.EQ.B2) THEN
               IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                  WRITE (REC,FMT=99988) ID(L), K, J, J, B1, BIGBND
                  CALL X04BAY(NERR,4,REC)
               END IF
            ELSE
               IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                  WRITE (REC,FMT=99987) ID(L), K, J, B1, J, B2
                  CALL X04BAY(NERR,3,REC)
               END IF
            END IF
         END IF
   40 CONTINUE
      IF (NERROR.GT.0) THEN
         IERR = 6
         GO TO 120
      END IF
C
      IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99985)
         CALL X04BAY(NOUT,2,REC)
         CALL A00AAF
         CALL H02BBT(M,N,TOLFES,BIGBND,TOLIV,MAXNOD,INTFST,MAXDPT,
     *               MSGLVL,ITMAX)
         IF (LRWORK.EQ.1 .AND. LIWORK.EQ.1) THEN
            WRITE (REC,FMT=99996) KLWORK, KRWORK, KIWORK
            IERR = 6
            GO TO 120
         ELSE
            WRITE (REC,FMT=99995) MAXDPT, LRWORK, LIWORK
            CALL X04BAY(NOUT,2,REC)
            WRITE (REC,FMT=99996) KLWORK, KRWORK, KIWORK
            CALL X04BAY(NOUT,2,REC)
         END IF
      END IF
      IF (ITMAX1.EQ.-11111) MSGLVL = MSG
      DO 60 I = 1, N
         IF (INTVAR(I).GT.0) GO TO 80
   60 CONTINUE
      IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
         WRITE (REC,FMT=99989)
         CALL X04BAY(NERR,2,REC)
      END IF
      IERR = 6
      GO TO 120
   80 CONTINUE
C
C     INITIALIZE REAL AND INTEGER WORKSPACE TO ZERO
C
      CALL F06FBF(LRWORK,RZERO,RWORK,1)
      CALL F06DBF(LIWORK,IZERO,IWORK,1)
C
C     CALCULATE NODEL AND NONOD
C
      NON1 = MAXDPT
      NON2 = (LRWORK-2*MIN(N,NCLIN+1)**2-13*N-12*NCLIN)/(N+2)
      IF (NON2.LE.0) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99997) LRWORK
            CALL X04BAY(NERR,2,REC)
         END IF
         IERR = 6
         GO TO 120
      END IF
      NONOD = MIN(NON1,NON2)
C
      IF (NONOD.LT.MAXDPT) THEN
         WRITE (REC,FMT=99997) LRWORK
         IERR = 8
         GO TO 120
      END IF
  100 CONTINUE
      NODEL = (LIWORK-3*N-(5+NCTOTL)*NONOD-NCTOTL-4)/5
C
      IF (NODEL.LT.4*NONOD) THEN
         IF (NONOD.GT.3*MAXDPT) THEN
            NONOD = NONOD/2
            GO TO 100
         ELSE
            IERR = 8
            WRITE (REC,FMT=99997) LIWORK
            GO TO 120
         END IF
      END IF
C
C     WORKSPACE REQUIREMENTS OF E04MFF
C
      LIWRK = 2*N + 3
      LRWRK = 2*MIN(N,NCLIN+1)**2 + 7*N + 5*NCLIN
C
      IOPTCL = NONOD + 1
C
C     CALCULATE WORK AREA OFFSETS
C
      BVINDX = LIWRK + 1
      SPLTVL = BVINDX + NODEL
      FATHER = SPLTVL + NODEL
      LFTSON = FATHER + NODEL
      RITSON = LFTSON + NODEL
      NXANOD = RITSON + NODEL
      ACTNOD = NXANOD + NONOD + 1
      STACK = ACTNOD + NONOD
      PARNOD = STACK + NONOD
      SOLVED = PARNOD + NONOD
      ISTATE = SOLVED + NONOD
C
      BLO = LRWRK + 1
      BUP = BLO + NCTOTL
      CLAMDA = BUP + NCTOTL
      BL1 = CLAMDA + NCTOTL
      BU1 = BL1 + NCTOTL
      CLAM = BU1 + NCTOTL
      AOPTVL = CLAM + NCTOTL
      FEASPT = AOPTVL + NONOD
      SOLN = FEASPT + NONOD
      AX = SOLN + N*NONOD
C
      CALL H02BBZ(ITMAX,MSGLVL,MAXNOD,N,NCLIN,NCTOTL,A,LDA,BL,BU,INTVAR,
     *            CVEC,TOLIV,TOLFES,BIGBND,X,OBJMIP,IWORK(BVINDX),
     *            IWORK(SPLTVL),IWORK(FATHER),IWORK(LFTSON),
     *            IWORK(RITSON),IWORK(PARNOD),IWORK(NXANOD),
     *            IWORK(ACTNOD),IWORK(STACK),IWORK(SOLVED),IWORK(ISTATE)
     *            ,RWORK(BLO),RWORK(BUP),RWORK(CLAMDA),RWORK(BL1),
     *            RWORK(BU1),RWORK(CLAM),RWORK(AOPTVL),RWORK(FEASPT),
     *            RWORK(SOLN),RWORK(AX),IWORK,LIWRK,RWORK,LRWRK,NONOD,
     *            IOPTCL,NODEL,INTFST,NODKNT,MAXDPT,IERR,REC)
C
      IF (IERR.EQ.0 .OR. IERR.EQ.7 .OR. IERR.EQ.9) THEN
C
C        SAVE THE INTEGER SOLUTION
C
         CALL H02BBP(RWORK(BL1),RWORK(BU1),RWORK(CLAM),IWORK(ISTATE),
     *               NCTOTL,IOPTCL,IWORK,LIWORK,RWORK,LRWORK)
      END IF
C
  120 IF (ITMAX1.EQ.-11111) MSGLVL = MSGLV1
      IF (IERR.EQ.6) WRITE (REC,FMT=99986)
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,2,REC)
      RETURN
C
99999 FORMAT (' ** On entry, LDA.lt.max(1,M):',/'    LDA =',I16,' M =',
     *       I16)
99998 FORMAT (' ** On entry, element',I16,' of INTVAR is not equal to ',
     *       '0 or 1',/'    INTVAR(',I16,') =',I16)
99997 FORMAT (' ** MAXDPT is too small to solve the problem:',/'    MA',
     *       'XDPT =',I16)
99996 FORMAT (' ** Workspace required with MAXDPT =',I6,': LRWORK =',I6,
     *       '  LIWORK =',I6,/)
99995 FORMAT (/' ** Workspace provided with MAXDPT =',I6,': LRWORK =',
     *       I6,'  LIWORK =',I6)
99994 FORMAT (' ** On entry, the dimension of RWORK is too small: LRWO',
     *       'RK =',I16,'.',/'    RWORK must be of dimension (at least)'
     *       ,I16)
99993 FORMAT (' ** On entry, the dimension of IWORK is too small: LIWO',
     *       'RK =',I16,'.',/'    IWORK must be of dimension (at least)'
     *       ,I16)
99992 FORMAT (' ** On entry, N.le.0:',/'    N =',I16)
99991 FORMAT (' ** On entry, M.lt.0:',/'    M =',I16)
99990 FORMAT (' ** On entry, MAXDPT.lt.2:',/'    MAXDPT =',I16)
99989 FORMAT (' ** On entry, all elements of INTVAR are set to zero',
     *       /'    i.e. there are no integer variables in the problem')
99988 FORMAT (/' ** On entry, the equal bounds on  ',A5,I3,'  are infi',
     *       'nite (because',/'    BL(',I4,').eq.beta and BU(',I4,').e',
     *       'q.beta, but abs(beta).ge.bigbnd):',/'    beta =',G16.7,
     *       ' bigbnd =',G16.7)
99987 FORMAT (/' ** On entry, the bounds on  ',A5,I3,'  are inconsiste',
     *       'nt:',/'    BL(',I4,') =',G16.7,'   BU(',I4,') =',G16.7)
99986 FORMAT (/' ** An input parameter is invalid.  Problem abandoned.')
99985 FORMAT (/' *** H02BBF')
      END
