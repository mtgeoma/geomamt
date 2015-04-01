      SUBROUTINE H02BBZ(ITMAX,MSGLVL,MAXNOD,N,NCLIN,NCTOTL,A,LDA,BL,BU,
     *                  INTVAR,CVEC,TOLIV,TOLFES,BIGBND,X,OBJMIP,BVINDX,
     *                  SPLTVL,FATHER,LFTSON,RITSON,PARNOD,NXANOD,
     *                  ACTNOD,STACK,SOLVED,ISTATE,BLO,BUP,CLAMDA,BL1,
     *                  BU1,CLAM,AOPTVL,FEASPT,SOLN,AX,IWORK,LIWORK,
     *                  RWORK,LRWORK,NONOD,IOPTCL,NODEL,INTFST,NODKNT,
     *                  MAXDPT,IERR,REC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-929 (APR 1991).
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1992.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, OBJMIP, TOLFES, TOLIV
      INTEGER           IERR, INTFST, IOPTCL, ITMAX, LDA, LIWORK,
     *                  LRWORK, MAXDPT, MAXNOD, MSGLVL, N, NCLIN,
     *                  NCTOTL, NODEL, NODKNT, NONOD
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AOPTVL(NONOD), AX(*), BL(NCTOTL),
     *                  BL1(NCTOTL), BLO(NCTOTL), BU(NCTOTL),
     *                  BU1(NCTOTL), BUP(NCTOTL), CLAM(NCTOTL),
     *                  CLAMDA(NCTOTL), CVEC(N), FEASPT(NONOD),
     *                  RWORK(LRWORK), SOLN(N,NONOD), X(N)
      INTEGER           ACTNOD(NONOD), BVINDX(NODEL), FATHER(NODEL),
     *                  INTVAR(N), ISTATE(NCTOTL,IOPTCL), IWORK(LIWORK),
     *                  LFTSON(NODEL), NXANOD(0:NONOD), PARNOD(NONOD),
     *                  RITSON(NODEL), SOLVED(NONOD), SPLTVL(NODEL),
     *                  STACK(NONOD)
      CHARACTER*80      REC(4)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  DIFF, FINITY, OBJLP, RHSVAL, RLOVAL, RUPVAL,
     *                  ZLBAR, ZUBAR
      INTEGER           FREANO, FREDEL, I, IDEPTH, IFAIL, IHED, II,
     *                  IMAX, IPAREN, ISTAT, ITER, ITLIM, IVAR, J, JJ,
     *                  JM1, K, KM1, L, M, MSGL1, NODTST, NT, Q
      CHARACTER*16      STR
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      LOGICAL           H02BBU
      EXTERNAL          X02AMF, H02BBU
C     .. External Subroutines ..
      EXTERNAL          E04MFF, E04MHF, H02BBV, H02BBW, H02BBX, H02BBY,
     *                  X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX, NINT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      FINITY = 1.0D+0/X02AMF()
      CALL E04MHF('Nolist')
      II = 0
      WRITE (STR,FMT='(I10)') II
      CALL E04MHF('Print Level = '//STR)
C
C     INITIALIZE CONSTRUCTS
C
      WRITE (STR,FMT='(D16.8)') TOLFES
      CALL E04MHF('Feasibility Tolerance = '//STR)
      WRITE (STR,FMT='(D16.8)') BIGBND
      CALL E04MHF('Infinite Bound Size = '//STR)
C
C     INITIALIZE FREE DELTA NODE LIST
C
      DO 20 K = 2, NODEL - 1
         RITSON(K) = K + 1
   20 CONTINUE
      FREDEL = 2
      RITSON(NODEL) = 0
C
      IPAREN = 1
      ISTAT = 0
      IHED = 0
      IDEPTH = 0
      ITLIM = 0
      OBJMIP = 0.0D+0
C
C     INITIALIZE FREE ACTIVE NODE LIST
C
      DO 40 K = 2, NONOD - 1
         NXANOD(K) = K + 1
   40 CONTINUE
      FREANO = 2
      NXANOD(NONOD) = 0
      NODKNT = 0
      MSGL1 = MSGLVL
      NODTST = 1
C
C     SET LOWER & UPPER BOUNDS ON CONSTRAINTS
C
      DO 60 K = N + 1, NCTOTL
         BLO(K) = BL(K)
         BUP(K) = BU(K)
   60 CONTINUE
C
C     SET RELAXED LP PROBLEM AS FIRST ACTIVE NODE
C
      K = 1
      NXANOD(0) = 1
      NXANOD(1) = 0
      ACTNOD(1) = 1
      LFTSON(1) = 0
      RITSON(1) = 0
      ZUBAR = FINITY
C     COLD = .TRUE.
      CALL E04MHF('Cold Start')
C
C     FIND A LOWER BOUND - PREPARE NODE AS AN LP
C
C     FIND THE PATH TO THE ROOT
C
   80 CONTINUE
      I = 0
      L = ACTNOD(K)
C
  100 CONTINUE
      IF (L.EQ.1) GO TO 120
      I = I + 1
      STACK(I) = L
      L = FATHER(L)
      GO TO 100
  120 CONTINUE
C
      IDEPTH = I
C
C     COPY ORIGINAL UPPER & LOWER VARIABLE BOUNDS
C
      DO 140 JJ = 1, N
         BLO(JJ) = BL(JJ)
         BUP(JJ) = BU(JJ)
  140 CONTINUE
C
C     APPLY DELTAS WHILE TRAVERSING THE TREE
C
      L = 1
      DO 160 JJ = I, 1, -1
         M = STACK(JJ)
         IF (LFTSON(L).EQ.M) THEN
            BUP(BVINDX(L)) = SPLTVL(L)
         ELSE
            BLO(BVINDX(L)) = SPLTVL(L) + 1
         END IF
         L = M
  160 CONTINUE
C
      NODKNT = NODKNT + 1
      IF (NODKNT.LT.2) GO TO 180
      RUPVAL = BUP(IVAR)
      RLOVAL = BLO(IVAR)
      RHSVAL = X(IVAR)
C
C     CALL THE LP SOLVER
C
  180 IFAIL = 1
      CALL E04MFF(N,NCLIN,A,LDA,BLO,BUP,CVEC,ISTATE(1,K),X,ITER,OBJLP,
     *            AX,CLAMDA,IWORK,LIWORK,RWORK,LRWORK,IFAIL)
C
      IF (MSGL1.GT.1 .AND. IFAIL.LT.2) THEN
         WRITE (REC,FMT=99996) OBJLP
         CALL X04BAY(NOUT,2,REC)
      END IF
      IF (MSGL1.EQ.10 .OR. MSGL1.EQ.20) THEN
         II = 5
         CALL H02BBW(N,NCTOTL,A,LDA,BLO,BUP,X,OBJMIP,K,ISTATE,CLAMDA,
     *               NODKNT,II,BIGBND)
      END IF
      IF (MSGL1.GT.1 .AND. IFAIL.LT.2) THEN
         WRITE (REC,FMT=99995)
         CALL X04BAY(NOUT,4,REC)
      END IF
      MSGL1 = 0
      IF (NODKNT.LT.2) THEN
         IF (IFAIL.EQ.2) THEN
            WRITE (REC,FMT=99988)
            IERR = 2
            RETURN
         END IF
         IF (IFAIL.EQ.3) THEN
            WRITE (REC,FMT=99987)
            IERR = 3
            RETURN
         END IF
         IF (IFAIL.EQ.4) THEN
            WRITE (REC,FMT=99986)
            IERR = 4
            RETURN
         END IF
      END IF
      IF (IFAIL.EQ.3 .OR. IFAIL.EQ.6) THEN
         OBJLP = FINITY
         IF (MSGLVL.GT.1) THEN
            IF (RLOVAL.LE.-BIGBND) THEN
               WRITE (REC,FMT=99984) NODKNT, IPAREN, IVAR, RHSVAL,
     *           RUPVAL, X(IVAR), IDEPTH
            ELSE IF (RUPVAL.GE.BIGBND) THEN
               WRITE (REC,FMT=99983) NODKNT, IPAREN, IVAR, RHSVAL,
     *           RLOVAL, X(IVAR), IDEPTH
            ELSE IF (RLOVAL.LE.-BIGBND .AND. RUPVAL.GE.BIGBND) THEN
               WRITE (REC,FMT=99982) NODKNT, IPAREN, IVAR, RHSVAL,
     *           X(IVAR), IDEPTH
            ELSE
               WRITE (REC,FMT=99985) NODKNT, IPAREN, IVAR, RHSVAL,
     *           RLOVAL, RUPVAL, X(IVAR), IDEPTH
            END IF
            CALL X04BAF(NOUT,REC(1))
         END IF
C
      ELSE IF (IFAIL.GE.2) THEN
         IF (IFAIL.EQ.4 .AND. ITLIM.EQ.0) THEN
            ITLIM = 1
            IMAX = MAX(100,10*(N+NCLIN))
            WRITE (STR,FMT='(I10)') IMAX
            CALL E04MHF('Iteration Limit = '//STR)
            CALL E04MHF('Warm Start')
            IF (MSGLVL.GT.1) THEN
               WRITE (REC,FMT=99999)
               CALL X04BAY(NOUT,2,REC)
            END IF
            GO TO 180
         END IF
         OBJLP = FINITY
         IF (MSGLVL.GT.1) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAY(NOUT,2,REC)
         END IF
         IF (IFAIL.EQ.4) WRITE (REC,FMT=99997)
         IERR = 7
         ISTAT = 2
      END IF
      AOPTVL(K) = OBJLP
      IF (NODKNT.LT.2) THEN
         IF (MSGLVL.GT.1) CALL H02BBY
      ELSE
         IF (IFAIL.LT.2 .AND. ZUBAR.GT.AOPTVL(K)) THEN
            IF (MSGLVL.GT.1) THEN
               IF (RLOVAL.LE.-BIGBND) THEN
                  WRITE (REC,FMT=99980) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, RUPVAL, X(IVAR), IDEPTH
               ELSE IF (RUPVAL.GE.BIGBND) THEN
                  WRITE (REC,FMT=99979) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, RLOVAL, X(IVAR), IDEPTH
               ELSE IF (RLOVAL.LE.-BIGBND .AND. RUPVAL.GE.BIGBND) THEN
                  WRITE (REC,FMT=99978) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, X(IVAR), IDEPTH
               ELSE
                  WRITE (REC,FMT=99981) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, RLOVAL, RUPVAL, X(IVAR), IDEPTH
               END IF
               CALL X04BAF(NOUT,REC(1))
            END IF
         END IF
      END IF
C
C     OBTAIN Z-LOWER-BAR
C
      ZLBAR = FINITY
      I = 0
  200 CONTINUE
      I = NXANOD(I)
      IF (I.EQ.0) GO TO 220
      IF (AOPTVL(I).LT.ZLBAR) ZLBAR = AOPTVL(I)
      GO TO 200
  220 CONTINUE
C
C     CHECK WHETHER SOLN IS WITHIN BOUNDS; WHETHER FEASIBLE
C
      IF (ZUBAR.LE.AOPTVL(K)) THEN
         IF (IFAIL.LT.2) THEN
            IF (MSGLVL.GT.1) THEN
               IF (RLOVAL.LE.-BIGBND) THEN
                  WRITE (REC,FMT=99976) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, RUPVAL, X(IVAR), IDEPTH
               ELSE IF (RUPVAL.GE.BIGBND) THEN
                  WRITE (REC,FMT=99975) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, RLOVAL, X(IVAR), IDEPTH
               ELSE IF (RLOVAL.LE.-BIGBND .AND. RUPVAL.GE.BIGBND) THEN
                  WRITE (REC,FMT=99974) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, X(IVAR), IDEPTH
               ELSE
                  WRITE (REC,FMT=99977) NODKNT, IPAREN, OBJLP, IVAR,
     *              RHSVAL, RLOVAL, RUPVAL, X(IVAR), IDEPTH
               END IF
               CALL X04BAF(NOUT,REC(1))
            END IF
         END IF
C
C        REMOVE ACTIVE NODE
C
         CALL H02BBX(KM1,K,NXANOD,FREANO,ACTNOD,FATHER,RITSON,LFTSON,
     *               FREDEL,NODTST)
      ELSE IF (H02BBU(X,N,INTVAR,TOLIV)) THEN
         ZUBAR = AOPTVL(K)
         IF (MSGLVL.GT.1) THEN
            WRITE (REC,FMT=99994)
            CALL X04BAY(NOUT,2,REC)
         END IF
C
         IHED = 1
         DO 240 II = 1, N
            FEASPT(II) = X(II)
  240    CONTINUE
         DO 260 II = 1, NCTOTL
            ISTATE(II,IOPTCL) = ISTATE(II,K)
            BL1(II) = BLO(II)
            BU1(II) = BUP(II)
            CLAM(II) = CLAMDA(II)
  260    CONTINUE
C
C        STOP IF FIRST INTEGER REQUIRED
C
         IF (INTFST.GT.0) THEN
            OBJMIP = ZUBAR
            ISTAT = 1
            IF (MSGLVL.GE.15) NODKNT = -NODKNT
            IF (MSGLVL.GT.0) CALL H02BBW(N,NCTOTL,A,LDA,BL1,BU1,X,
     *                                   OBJMIP,IOPTCL,ISTATE,CLAM,
     *                                   NODKNT,ISTAT,BIGBND)
            IERR = 0
            RETURN
         END IF
C
C        SOLN WAS FEASIBLE - REMOVE OBSOLETE ACTIVE NODES
C
         I = 0
  280    CONTINUE
         J = I
         I = NXANOD(J)
         IF (I.EQ.0) GO TO 300
         IF (ZUBAR.LE.AOPTVL(I)) THEN
C
C           REMOVE OBSOLETE ACTIVE NODES
C
            CALL H02BBX(J,I,NXANOD,FREANO,ACTNOD,FATHER,RITSON,LFTSON,
     *                  FREDEL,NODTST)
            I = J
         END IF
         GO TO 280
  300    CONTINUE
      ELSE
C
C        NODE WILL REMAIN ACTIVE
C
         IF (IDEPTH.GE.MAXDPT) THEN
            IERR = 8
            WRITE (REC,FMT=99972)
            RETURN
         END IF
         DO 320 II = 1, N
            SOLN(II,K) = X(II)
  320    CONTINUE
         SOLVED(K) = 1
C
         PARNOD(K) = NODKNT
C
      END IF
C
C     CHECK WHETHER ACTIVE NODE LIST IS EMPTY
C
      IF (NXANOD(0).EQ.0) THEN
         IF (ZUBAR.LT.FINITY) THEN
            OBJMIP = ZUBAR
            DO 340 II = 1, N
               X(II) = FEASPT(II)
  340       CONTINUE
C
            IF (MSGLVL.GT.1) THEN
               WRITE (REC,FMT=99993)
               CALL X04BAF(NOUT,REC(1))
            END IF
C
            IF (MSGLVL.GE.15) NODKNT = -NODKNT
            IF (MSGLVL.GT.0) CALL H02BBW(N,NCTOTL,A,LDA,BL1,BU1,X,
     *                                   OBJMIP,IOPTCL,ISTATE,CLAM,
     *                                   NODKNT,ISTAT,BIGBND)
C
         ELSE
            OBJMIP = FINITY
            IERR = 1
            WRITE (REC,FMT=99992)
            RETURN
         END IF
         IERR = 0
         RETURN
      END IF
C
C     SELECT A NEW ACTIVE NODE TO WORK ON
C
      IF (IHED.GT.0 .AND. MSGLVL.GT.1) THEN
         CALL H02BBY
         IHED = 0
      END IF
C
C     SELECT A NEW PROBLEM
C
      IF (IDEPTH.GT.MAXDPT) THEN
         IERR = 8
         WRITE (REC,FMT=99972)
         RETURN
      END IF
      CALL H02BBV(NXANOD,JM1,J,ZLBAR,AOPTVL,NODTST,NONOD)
      DO 360 II = 1, N
         X(II) = SOLN(II,J)
  360 CONTINUE
C
C     CHECK WHETHER ACTIVE NODE HAS BEEN SOLVED BUT NOT BRANCHED
C
      IF (SOLVED(J).EQ.1) THEN
         IPAREN = PARNOD(J)
         I = 1
  380    CONTINUE
         DIFF = ABS(SOLN(I,J)-NINT(SOLN(I,J)))
         IF (INTVAR(I).EQ.1 .AND. DIFF.GT.TOLIV) GO TO 400
         I = I + 1
         GO TO 380
  400    CONTINUE
C
         IVAR = I
C
         L = ACTNOD(J)
         BVINDX(L) = I
         SPLTVL(L) = INT(SOLN(I,J))
         IF (SOLN(I,J).LT.0.0D+0) SPLTVL(L) = SPLTVL(L) - 1
C
C        BRANCH - ADD LEFT SON
C
         IF (FREDEL.EQ.0) THEN
            IERR = 8
            WRITE (REC,FMT=99973)
            RETURN
         END IF
         M = FREDEL
         FREDEL = RITSON(M)
         LFTSON(L) = M
         FATHER(M) = L
         LFTSON(M) = 0
         RITSON(M) = 0
         IF (FREANO.EQ.0) THEN
            IERR = 8
            WRITE (REC,FMT=99973)
            RETURN
         END IF
C
C        TEST FOR MAXIMUM NODE COUNT
C
         IF (MAXNOD.GT.0 .AND. MAXNOD.LE.NODKNT+1) THEN
            IF (MSGLVL.GT.0) THEN
               WRITE (REC,FMT=99991)
               CALL X04BAY(NOUT,4,REC)
            END IF
            IF (ZUBAR.LT.FINITY) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99990)
                  CALL X04BAY(NOUT,2,REC)
               END IF
               DO 420 II = 1, N
                  X(II) = FEASPT(II)
  420          CONTINUE
               OBJMIP = ZUBAR
               ISTAT = 2
               IF (MSGLVL.GE.15) NODKNT = -NODKNT
               IF (MSGLVL.GT.0) THEN
                  NT = NODKNT + 1
                  CALL H02BBW(N,NCTOTL,A,LDA,BL1,BU1,X,OBJMIP,IOPTCL,
     *                        ISTATE,CLAM,NT,ISTAT,BIGBND)
               END IF
               WRITE (REC,FMT=99971)
               IERR = 9
               RETURN
            ELSE
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99989)
                  CALL X04BAY(NOUT,2,REC)
               END IF
               OBJMIP = FINITY
               WRITE (REC,FMT=99971)
               IERR = 10
               RETURN
            END IF
         END IF
C
         NODTST = NODTST + 1
         Q = FREANO
         FREANO = NXANOD(Q)
         NXANOD(Q) = NXANOD(0)
         NXANOD(0) = Q
         IF (JM1.EQ.0) JM1 = Q
         ACTNOD(Q) = M
C
         AOPTVL(Q) = AOPTVL(J)
         DO 440 II = 1, N
            SOLN(II,Q) = SOLN(II,J)
  440    CONTINUE
C
         SOLVED(Q) = 0
         DO 460 II = 1, NCTOTL
            ISTATE(II,Q) = ISTATE(II,J)
  460    CONTINUE
         ISTATE(I,Q) = -1
C
C        ADD RIGHT SON & DE-ACTIVATE FATHER NODE
C
         IF (FREDEL.EQ.0) THEN
            IERR = 8
            WRITE (REC,FMT=99973)
            RETURN
         END IF
         M = FREDEL
         FREDEL = RITSON(M)
         RITSON(L) = M
         FATHER(M) = L
         LFTSON(M) = 0
         RITSON(M) = 0
         ACTNOD(J) = M
C
         SOLVED(J) = 0
         ISTATE(I,J) = -2
C
      END IF
      KM1 = JM1
      K = J
      GO TO 80
C
99999 FORMAT (/' ** Doubling iteration limit')
99998 FORMAT (/' ** Cannot find optimum LP, terminate the search of th',
     *       'is branch')
99997 FORMAT (/' ** Search of a branch was terminated due to iteration',
     *       ' limit')
99996 FORMAT (/'   *** Optimum LP solution ***',G16.7)
99995 FORMAT (//'   *** Start of tree search ***',/)
99994 FORMAT ('      *** Integer solution ***',/)
99993 FORMAT ('     *** End of tree search ***')
99992 FORMAT (/' ** The problem does not have a feasible integer solut',
     *       'ion.')
99991 FORMAT (//' ** Number of nodes exceeds the maximum **',/)
99990 FORMAT (' ** The best integer solution **',/)
99989 FORMAT (' ** No integer solution so far **',/)
99988 FORMAT (/' ** The LP solution is unbounded.')
99987 FORMAT (/' ** The LP does not have a feasible solution.')
99986 FORMAT (/' ** Iteration limit reached without finding a solution.'
     *       )
99985 FORMAT (1X,I4,2X,I4,4X,'No Feas Soln',2X,I3,3X,G9.3,1X,G9.3,2X,
     *       G9.3,1X,G9.3,1X,I3)
99984 FORMAT (1X,I4,2X,I4,4X,'No Feas Soln',2X,I3,3X,G9.3,2X,'None',6X,
     *       G9.3,1X,G9.3,1X,I3)
99983 FORMAT (1X,I4,2X,I4,4X,'No Feas Soln',2X,I3,3X,G9.3,1X,G9.3,3X,
     *       'None',5X,G9.3,1X,I3)
99982 FORMAT (1X,I4,2X,I4,4X,'No Feas Soln',2X,I3,3X,G9.3,2X,'None',8X,
     *       'None',5X,G9.3,1X,I3)
99981 FORMAT (1X,I4,2X,I4,4X,G10.3,4X,I3,3X,G9.3,1X,G9.3,2X,G9.3,1X,
     *       G9.3,1X,I3)
99980 FORMAT (1X,I4,2X,I4,4X,G10.3,4X,I3,3X,G9.3,2X,'None',6X,G9.3,1X,
     *       G9.3,1X,I3)
99979 FORMAT (1X,I4,2X,I4,4X,G10.3,4X,I3,3X,G9.3,1X,G9.3,3X,'None',5X,
     *       G9.3,1X,I3)
99978 FORMAT (1X,I4,2X,I4,4X,G10.3,4X,I3,3X,G9.3,2X,'None',7X,'None',
     *       5G9.3,1X,I3)
99977 FORMAT (1X,I4,2X,I4,4X,G10.3,' CO',1X,I3,3X,G9.3,1X,G9.3,2X,G9.3,
     *       1X,G9.3,1X,I3)
99976 FORMAT (1X,I4,2X,I4,4X,G10.3,' CO',1X,I3,3X,G9.3,2X,'None',6X,
     *       G9.3,1X,G9.3,1X,I3)
99975 FORMAT (1X,I4,2X,I4,4X,G10.3,' CO',1X,I3,3X,G9.3,1X,G9.3,3X,
     *       'None',5X,G9.3,1X,I3)
99974 FORMAT (1X,I4,2X,I4,4X,G10.3,' CO',1X,I3,3X,G9.3,2X,'None',7X,
     *       'None',5X,G9.3,1X,I3)
99973 FORMAT (/' ** Not enough workspace to solve problem.')
99972 FORMAT (/' ** Maximum depth value is too small. Increase MAXDPT ',
     *       'and rerun H02BBF.')
99971 FORMAT (/' ')
      END
