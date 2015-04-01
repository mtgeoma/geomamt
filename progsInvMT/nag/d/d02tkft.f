      SUBROUTINE D02TKF(FFUN,FJAC,GAFUN,GBFUN,GAJAC,GBJAC,GUESS,WORK,
     *                  IWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02TKF')
      INTEGER           FLAG, NUMBER, MLIMIT, ALTER
      PARAMETER         (FLAG=1,NUMBER=2,MLIMIT=3,ALTER=4)
C     .. Scalar Arguments ..
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*)
      INTEGER           IWORK(*)
C     .. Subroutine Arguments ..
      EXTERNAL          FFUN, FJAC, GAFUN, GAJAC, GBFUN, GBJAC, GUESS
C     .. Local Scalars ..
      DOUBLE PRECISION  ALEFT, ARIGHT, ERMX
      INTEGER           I, ICARE, IER, IERMX, IFLAG, IGUESS, IJERMX,
     *                  IOUT, IPRINT, IREAD, J, K, KD, LACCUM, LCOEF,
     *                  LDELDZ, LDELZ, LDF, LDGR, LDGZ, LDMZ, LDQDMZ,
     *                  LDQZ, LDSCL, LEREST, LERMX, LERR, LF, LFIXP, LG,
     *                  LINTEG, LIPMSH, LJTOL, LLDGR, LLTOL, LM, LMSHD,
     *                  LNEWXI, LPVTG, LPVTW, LRHS, LROOT, LSCL, LSDMZ,
     *                  LSLOPE, LSXI, LSZ, LTOL, LV, LVALST, LW, LWTERR,
     *                  LWTMSH, LXI, LXIOLD, LZ, LZVAL, LZVAL1, MATCH,
     *                  MAXORD, MMAX, MODE, MSTAR, MXFIXP, N, NDMZ, NEQ,
     *                  NFIXP, NLBC, NMAX, NOLD, NONLIN, NP1, NREC,
     *                  NTOL, NZ
      LOGICAL           SVMESH
C     .. Local Arrays ..
      DOUBLE PRECISION  ACOL(28,7), ASAVE(28,4), B(28), DUMMY(1), RHO(7)
      INTEGER           MSHINF(4)
      CHARACTER*80      REC(5)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02TKL, D02TKW, D02TKZ, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
      IER = 0
      NREC = 0
      IF (IWORK(1).NE.1) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,a/a)')
     *     ' ** It appears that the setup routine D02TVF has not been',
     *     ' called.', ' ** No solution will be computed.'
         GO TO 340
      END IF
C 
      IPRINT = IWORK(17)
      CALL X04ABF(0,IOUT)
C 
C  In case incorrect input data is detected, the program returns
C  immediately with IFLAG=-3.
C 
      NEQ = IWORK(2)
      IFLAG = -3
C 
C  Rename some of the parameters and set default values.
C 
      NONLIN = 1
      K = IWORK(3)
      N = IWORK(4)
      IF (N.EQ.0) N = 5
C  IREAD = IPAR(8)
      IREAD = 1
C  IGUESS = IPAR(9)
      IGUESS = IWORK(10)
      IF (NONLIN.EQ.0 .AND. IGUESS.EQ.1) IGUESS = 0
      IF (IGUESS.GE.2 .AND. IREAD.EQ.0) IREAD = 1
C  ICARE = IPAR(10)
      ICARE = 0
C  Problem sizes
      NTOL = NEQ
      NFIXP = IWORK(8)
      MXFIXP = IWORK(20)
      MSTAR = IWORK(7)
      MMAX = IWORK(6)
      MAXORD = MMAX
      IF (K.EQ.0) K = MAX(MMAX+1,5-MMAX)
      KD = K*NEQ
C 
C  Print the input data for checking.
C 
      IF (IPRINT.LT.1) THEN
         IF (NONLIN.GT.0) THEN
            WRITE (REC(1),FMT='(a,i4,a)')
     *        'DEB Solving nonlinear system of ', NEQ,
     *        ' equations with orders:'
         ELSE
            WRITE (REC(1),FMT='(a,i4,a)')
     *        'DEB Solving linear system of ', NEQ,
     *        ' equations with orders:'
         END IF
         CALL X04BAF(IOUT,REC(1))
         DO 20 I = 1, NEQ, 20
            WRITE (REC(1),FMT='(a,20i3)') 'DEB ',
     *        (IWORK(20+J),J=I,MIN(I+19,NEQ))
            CALL X04BAF(IOUT,REC(1))
   20    CONTINUE
         WRITE (REC(1),FMT='(a,i4,a)') 'DEB There are ', NFIXP,
     *     ' fixed internal mesh points:'
         CALL X04BAF(IOUT,REC(1))
         DO 40 I = 1, NFIXP, 5
            WRITE (REC(1),FMT='(a,5e13.5)') 'DEB ',
     *        (WORK(1+J),J=I,MIN(I+4,NFIXP))
            CALL X04BAF(IOUT,REC(1))
   40    CONTINUE
         WRITE (REC(1),FMT='(a,i1,a)') 'DEB There are ', K,
     *     ' collocation points per subinterval'
         CALL X04BAF(IOUT,REC(1))
         WRITE (REC(1),FMT='(a)')
     *     'DEB Components of Z requiring tolerances are:'
         CALL X04BAF(IOUT,REC(1))
         DO 60 I = 1, NEQ, 20
            WRITE (REC(1),FMT='(a,20i3)') 'DEB ',
     *        (IWORK(20+2*NEQ+J),J=I,MIN(I+19,NEQ))
            CALL X04BAF(IOUT,REC(1))
   60    CONTINUE
         WRITE (REC(1),FMT='(a)') 'DEB Corresponding tolerances are :'
         CALL X04BAF(IOUT,REC(1))
         DO 80 I = 1, NEQ, 5
            WRITE (REC(1),FMT='(a,5e13.5)') 'DEB ',
     *        (WORK(1+MXFIXP+J),J=I,MIN(I+4,NEQ))
            CALL X04BAF(IOUT,REC(1))
   80    CONTINUE
         IF (IGUESS.GE.2) THEN
            WRITE (REC(1),FMT='(a)')
     *        'DEB Continutation mesh and solution are user supplied'
            CALL X04BAF(IOUT,REC(1))
         END IF
         IF (IREAD.EQ.2) THEN
            WRITE (REC(1),FMT='(a)') 'DEB No adaptive mesh selection'
            CALL X04BAF(IOUT,REC(1))
         END IF
      END IF
C 
      NLBC = IWORK(9)
C 
C  Set limits on iterations and initialize counters.
C  limit = maximum number of Newton iterations per mesh.
C  See subroutine  newmsh  for the roles of MSHINF(MLIMIT), MSHINF(FLAG)
C  MSHINF(NUMBER), and MSHINF(ALTER) .
C 
      MSHINF(MLIMIT) = 3
      MSHINF(FLAG) = 0
      MSHINF(NUMBER) = 1
      MSHINF(ALTER) = 1
C 
      NMAX = IWORK(5)
C 
C  Generate pointers to break up WORK and IWORK.
C 
      LERMX = 1
      LFIXP = 2
      LTOL = LFIXP + MXFIXP
      LSXI = LTOL + NEQ
      LSZ = LSXI + NMAX + 1
      LSDMZ = LSZ + MSTAR*(NMAX+1)
      LCOEF = LSDMZ + KD*NMAX
      LNEWXI = LCOEF + 49
      LXI = LNEWXI + NMAX + 1
      LROOT = LXI + NMAX + 1
      LWTMSH = LROOT + NEQ
      LWTERR = LWTMSH + NEQ
      LDGZ = LWTERR + MSTAR
      LDF = LDGZ + MSTAR
      LF = LDF + NEQ*NEQ*MAXORD
      LZVAL = LF + NEQ
      LZVAL1 = LZVAL + MSTAR
      LMSHD = LZVAL1 + NEQ*MAXORD
      LERR = LMSHD + 2*NEQ
      LEREST = LERR + MSTAR
      LDGR = LEREST + MSTAR
      LLDGR = MAX(NLBC,MSTAR-NLBC)
      LG = LDGR + NEQ*MAXORD*LLDGR
      LXIOLD = LG + (2*NMAX+1)*MSTAR*MSTAR
      LW = LXIOLD + NMAX + 1
      LV = LW + KD**2*NMAX
      LZ = LV + MSTAR*KD*NMAX
      LDMZ = LZ + MSTAR*(NMAX+1)
      LDELZ = LDMZ + KD*NMAX
      LDELDZ = LDELZ + MSTAR*(NMAX+1)
      LDQZ = LDELDZ + KD*NMAX
      LDQDMZ = LDQZ + MSTAR*(NMAX+1)
      LRHS = LDQDMZ + KD*NMAX
      LVALST = LRHS + KD*NMAX + MSTAR
      LSLOPE = LVALST + 4*MSTAR*NMAX
      LACCUM = LSLOPE + NMAX
      LSCL = LACCUM + NMAX + 1
      LDSCL = LSCL + MSTAR*(NMAX+1)
C 
      LM = 21
      LJTOL = LM + NEQ
      LLTOL = LJTOL + NEQ
      LIPMSH = LLTOL + NEQ
      LPVTG = LIPMSH + NMAX + 1
      LPVTW = LPVTG + MSTAR*(NMAX+1)
      LINTEG = LPVTW + KD*NMAX
C 
C  If  IGUESS .GE. 2, move XIOLD, Z, and DMZ to their proper
C  locations in WORK.
C 
      IF (IGUESS.EQ.4) THEN
         ALEFT = WORK(LNEWXI)
         ARIGHT = WORK(LNEWXI+N)
      ELSE
         ALEFT = WORK(LXI)
         ARIGHT = WORK(LXI+N)
      END IF
      IF (IGUESS.GE.2) THEN
         NOLD = IWORK(11)
         NZ = MSTAR*(NOLD+1)
         NDMZ = KD*NOLD
         DO 100 I = 1, NZ
            WORK(LZ+I-1) = WORK(LSZ+I-1)
  100    CONTINUE
         DO 120 I = 1, NDMZ
            WORK(LDMZ+I-1) = WORK(LSDMZ+I-1)
  120    CONTINUE
         NP1 = NOLD + 1
         DO 140 I = 1, NP1
            WORK(LXIOLD+I-1) = WORK(LSXI+I-1)
  140    CONTINUE
         DO 160 I = 1, N + 1
            WORK(LXI+I-1) = WORK(LNEWXI+I-1)
  160    CONTINUE
      END IF
C 
C  Initialize collocation points, constants, mesh.
C 
      CALL D02TKZ(K,RHO,WORK(LCOEF),ACOL,B,ASAVE,IWORK(LM),NEQ,MAXORD,
     *            WORK(LWTERR),IWORK(LJTOL),IWORK(LLTOL),WORK(LROOT),
     *            WORK(LWTMSH),MSTAR,NTOL,WORK(LTOL))
      MODE = 3 + IREAD
      CALL D02TKL(MODE,ALEFT,ARIGHT,WORK(LXI),WORK(LXIOLD),N,NOLD,NMAX,
     *            DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,NFIXP,WORK(LFIXP),
     *            K,NEQ,IWORK(LM),MAXORD,MSTAR,WORK(LCOEF),IGUESS,ASAVE,
     *            NTOL,IWORK(LJTOL),IWORK(LLTOL),WORK(LWTMSH),
     *            WORK(LROOT),MSHINF,IOUT,IPRINT)
      NZ = MSTAR*(N+1)
      NDMZ = NEQ*K*N
C 
C  Determine first approximation, if the problem is nonlinear.
C 
      IF (IGUESS.GE.2) GO TO 240
      NP1 = N + 1
      DO 180 I = 1, NP1
         WORK(I+LXIOLD-1) = WORK(I+LXI-1)
  180 CONTINUE
      NOLD = N
      IF (NONLIN.EQ.0 .OR. IGUESS.EQ.1) GO TO 240
C 
C  System provides first approximation of the solution.
C  Choose Z(j) = 0  for j=1,...,mstar.
C 
      DO 200 I = 1, NZ
         WORK(LZ-1+I) = 0.D0
  200 CONTINUE
      DO 220 I = 1, NDMZ
         WORK(LDMZ-1+I) = 0.D0
  220 CONTINUE
C 
C  Ok! Let's do it!
C 
  240 CONTINUE
      ERMX = 0.0D0
      IERMX = 0
      IJERMX = 0
      IF (IGUESS.GE.2) IGUESS = 0
      CALL D02TKW(WORK(LXI),WORK(LXIOLD),WORK(LZ),WORK(LDMZ),WORK(LRHS),
     *            WORK(LDELZ),WORK(LDELDZ),WORK(LDQZ),WORK(LDQDMZ),
     *            WORK(LG),WORK(LW),WORK(LV),WORK(LDGZ),WORK(LDF),
     *            WORK(LF),WORK(LVALST),WORK(LSLOPE),WORK(LSCL),
     *            WORK(LDSCL),WORK(LACCUM),WORK(LMSHD),IWORK(LPVTG),
     *            IWORK(LINTEG),IWORK(LPVTW),NFIXP,WORK(LFIXP),IFLAG,
     *            FFUN,FJAC,GAFUN,GAJAC,GBFUN,GBJAC,GUESS,IWORK(LM),
     *            MAXORD,NEQ,MSTAR,MSHINF,WORK(LZVAL),WORK(LZVAL1),RHO,
     *            WORK(LCOEF),K,B,ACOL,ASAVE,ALEFT,ARIGHT,NLBC,KD,N,
     *            NOLD,NMAX,NZ,NDMZ,NONLIN,ICARE,IGUESS,NTOL,WORK(LTOL),
     *            WORK(LROOT),WORK(LWTMSH),WORK(LWTERR),IWORK(LJTOL),
     *            IWORK(LLTOL),WORK(LERR),WORK(LEREST),WORK(LDGR),LLDGR,
     *            ERMX,IERMX,IJERMX,WORK(LSXI),WORK(LSZ),WORK(LSDMZ),
     *            SVMESH,IOUT,IPRINT)
C 
C  Prepare output
C 
      WORK(LERMX) = ERMX
      IWORK(18) = IERMX
      IWORK(19) = IJERMX
      IWORK(4) = N
      IWORK(11) = N
C  Mesh change indicator
      IWORK(12) = 0
C 
      IF (IFLAG.EQ.1) THEN
         IER = 0
      ELSE IF (IFLAG.EQ.-1) THEN
         IER = 5
         NREC = 4
         WRITE (REC,FMT='(a/a,i6,a/a/a,a)')
     *   ' ** The expected number of subintervals required to continue '
     *     , ' ** the computation exceeds the maximum specified,',
     *     NMAX + 1, '.',
     *     ' ** Results have been generated which may be useful.',
     *     ' ** Try increasing this number or relaxing the error ',
     *     'requirements.'
      ELSE IF (IFLAG.EQ.0) THEN
         IER = 2
         NREC = 4
         WRITE (REC,FMT='(a/a/a/a)')
     *     ' ** Numerical singularity has been detected in the Jacobian'
     *     , ' ** used in the Newton iteration. No results have been',
     *     ' ** generated. Check the coding of the procedure arguments',
     *     ' ** FJAC, GAJAC and GBJAC.'
      ELSE IF (IFLAG.EQ.-2) THEN
         IF ( .NOT. SVMESH) THEN
            IER = 3
            NREC = 5
            WRITE (REC,FMT='(a/a/a/a/a)')
     *  ' ** All Newton iterations that have been attempted have failed'
     *        ,
     *     ' ** to converge. No results have been generated. Check the '
     *        ,
     *    ' ** coding of the procedure arguments FJAC, GAJAC and GBJAC.'
     *        ,
     *      ' ** Try to provide a better initial solution approximation'
     *        , ' ** in GUESS.'
         ELSE
            IER = 4
            NREC = 4
            WRITE (REC,FMT='(a/a/a/a)')
     *  ' ** A Newton iteration has failed to converge. The computation'
     *        ,
     *     ' ** has not succeeded but results have been returned for an'
     *        ,
     *  ' ** intermediate mesh on which convergence was achieved. These'
     *        , ' ** results should be treated with extreme caution.'
         END IF
      END IF
C 
C  Mesh point identifiers
C 
      IF (IER.EQ.0 .OR. IER.EQ.5 .OR. IER.EQ.4) THEN
         IWORK(LIPMSH) = 1
         DO 260 I = 1, N - 1
            IWORK(LIPMSH+I) = 2
  260    CONTINUE
         IWORK(LIPMSH+N) = 1
         DO 280 I = N + 1, NMAX + 1
            IWORK(LIPMSH+I) = -1
  280    CONTINUE
         MATCH = 0
         DO 300 I = 0, N
            IF (MATCH.LT.NFIXP) THEN
               IF (WORK(LFIXP+MATCH).EQ.WORK(LXI+I)) THEN
                  IWORK(LIPMSH+I) = 1
                  MATCH = MATCH + 1
               END IF
            END IF
  300    CONTINUE
         IF (IER.EQ.0 .OR. IER.EQ.5) THEN
            DO 320 I = 1, N, 2
               IWORK(LIPMSH+I) = 3
  320       CONTINUE
         END IF
      END IF
C 
  340 CONTINUE
      IWORK(16) = IER + 1
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C 
      RETURN
      END
