      SUBROUTINE E04NCZ(PRBTYP,NAMED,NAMES,LINOBJ,UNITQ,INFORM,ITER,
     *                  JINF,NCLIN,NCTOTL,NACTIV,NFREE,NRANK,NZ,NRZ,N,
     *                  LDA,LDR,ISTATE,KACTIV,KX,CTX,SSQ,SSQ1,SUMINF,
     *                  NUMINF,XNORM,BL,BU,A,CLAMDA,AX,FEATOL,R,X,W)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1076 (JUL 1993).
C     MARK 17 REVISED. IER-1586 (JUN 1995).
C
C     ******************************************************************
C     E04NCZ  is a subroutine for linearly constrained linear-least
C     squares.  On entry, it is assumed that an initial working set of
C     linear constraints and bounds is available.
C     The arrays ISTATE, KACTIV and KX will have been set accordingly
C     and the arrays T and ZY will contain the TQ factorization of
C     the matrix whose rows are the gradients of the active linear
C     constraints with the columns corresponding to the active bounds
C     removed.  the TQ factorization of the resulting (NACTIV by NFREE)
C     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
C     is reverse-triangular.
C
C     Values of ISTATE(J) for the linear constraints.......
C
C     ISTATE(J)
C     ---------
C          0    constraint J is not in the working set.
C          1    constraint J is in the working set at its lower bound.
C          2    constraint J is in the working set at its upper bound.
C          3    constraint J is in the working set as an equality.
C
C     Constraint J may be violated by as much as FEATOL(J).
C
C     Systems Optimization Laboratory, Stanford University.
C     This version of  E04NCZ  dated 14-Sep-1992.
C
C     Copyright  1984  Stanford University.
C
C     This material may be reproduced by or for the U.S. Government
C     pursuant to the copyright license under DAR clause 7-104.9(a)
C     (1979 Mar).
C
C     This material is based upon work partially supported by the
C     National Science Foundation under grants MCS-7926009 and
C     ECS-8012974; the Department of Energy Contract AM03-76SF00326, PA
C     No. DE-AT03-76ER72018; and the Army Research Office Contract DAA29
C     -79-C-0110.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      INTEGER           MSTALL, MREFN
      PARAMETER         (MSTALL=50,MREFN=1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CTX, SSQ, SSQ1, SUMINF, XNORM
      INTEGER           INFORM, ITER, JINF, LDA, LDR, N, NACTIV, NCLIN,
     *                  NCTOTL, NFREE, NRANK, NRZ, NUMINF, NZ
      LOGICAL           LINOBJ, NAMED, UNITQ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  CLAMDA(NCTOTL), FEATOL(NCTOTL), R(LDR,*), W(*),
     *                  X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, DTMAX,
     *                  DTMIN, EPSPT3, EPSPT5, EPSPT8, EPSPT9, TOLACT,
     *                  TOLFEA, TOLRNK
      INTEGER           IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LDT, LDZY, LENNAM, LFORMH, LINES1,
     *                  LINES2, LPROB, MSGLS, NCOLT, NN, NNCLIN, NOUT,
     *                  NPROB
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPSVLS(MXPARM), WMACH(15)
      INTEGER           IPADLS(19), IPSVLS(MXPARM), LOCLS(LENLS)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSRZZ, ALFA, ALFHIT, ATPHIT, BIGALF, CNORM,
     *                  CONDMX, CONDRZ, CONDT, CTP, DINKY, DRZMAX,
     *                  DRZMIN, ERR1, ERR2, FLMAX, GFNORM, GRZNRM,
     *                  GZNORM, OBJSIZ, PALFA, PNORM, RESNRM, ROWNRM,
     *                  TRULAM, WSSIZE
      INTEGER           IADD, IFIX, IREFN, IS, ISDEL, ITMAX, JADD,
     *                  JBIGST, JDEL, JMAX1, JSMLST, JTINY, KBIGST,
     *                  KDEL, KSMLST, LANORM, LAP, LCQ, LGQ, LHZ, LPX,
     *                  LRES, LRES0, LRLAM, LT, LWRK, LWTINF, LZY,
     *                  MSGLVL, NGQ, NPHASE, NRES, NSTALL, NVIOL
      LOGICAL           CONVRG, CYCLIN, ERROR, FIRSTV, HITCON, HITLOW,
     *                  NEEDFG, OVERFL, PRNT, ROWERR, SINGLR, STALL,
     *                  STATPT, UNBNDD, UNCON, UNITGZ, WEAK
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM)
      INTEGER           IPRMLS(MXPARM)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          E04MFK, E04NCH, E04NCJ, E04NCK, E04NCL, E04NCP,
     *                  E04NCQ, E04NCR, E04NCT, E04NCV, E04UCG, E04UPP,
     *                  F06FBF, F06FLF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (MSGLS,MSGLVL)
C     .. Save statement ..
      SAVE              /AX02ZA/, /DE04NC/, /EE04NC/
C     .. Executable Statements ..
C
C     Specify the machine-dependent parameters.
C
      FLMAX = WMACH(7)
C
      LANORM = LOCLS(2)
      LAP = LOCLS(3)
      LPX = LOCLS(4)
      LRES = LOCLS(5)
      LRES0 = LOCLS(6)
      LHZ = LOCLS(7)
      LGQ = LOCLS(8)
      LCQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LZY = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK = LOCLS(14)
C
C     Set up the adresses of the contiguous arrays  ( RES0, RES )
C     and  ( GQ, CQ ).
C
      NRES = 0
      IF (NRANK.GT.0) NRES = 2
      NGQ = 1
      IF (LINOBJ) NGQ = 2
C
C     Initialize.
C
      IREFN = 0
      ITER = 0
C
      IF (PRBTYP.EQ.'FP') THEN
         ITMAX = ITMAX2
      ELSE
         ITMAX = ITMAX1
      END IF
C
      JADD = 0
      JDEL = 0
      NPHASE = 1
      NSTALL = 0
      NUMINF = -1
      NRZ = 0
C
      ALFA = ZERO
      CONDMX = FLMAX
      DRZMAX = ONE
      DRZMIN = ONE
      SSQ = ZERO
C
      CYCLIN = .FALSE.
      ERROR = .FALSE.
      FIRSTV = .FALSE.
      PRNT = .TRUE.
      NEEDFG = .TRUE.
      STALL = .TRUE.
      UNCON = .FALSE.
      UNBNDD = .FALSE.
C
C     =================== start of the main loop =======================
C
C      cyclin = false
C      unbndd = false
C      error  = false
C      k      = 0
C
C      repeat
C            repeat
C                  compute Z'g,  print details of this iteration
C                  stat pt = (Z'g .eq. 0)
C                  if (not stat pt) then
C                     error =  k .ge. itmax
C                     if (not error) then
C                        compute p, alfa
C                        error = unbndd  or  cyclin
C                        if (not error) then
C                           k = k + 1
C                           x = x + alfa p
C                           if (feasible) update Z'g
C                           if necessary, add a constraint
C                        end if
C                     end if
C                  end if
C            until  stat pt  or  error
C
C            compute lam1, lam2, smllst
C            optmul =  smllst .gt. 0
C            if ( not (optmul .or. error) ) then
C                  delete an artificial or regular constraint
C            end if
C      until optmul  or  error
C
C     ==================================================================
C
C     REPEAT
C        REPEAT
   20 IF (NEEDFG) THEN
         IF (NRANK.GT.0) THEN
            RESNRM = DNRM2(NRANK,W(LRES),1)
            SSQ = HALF*(SSQ1**2+RESNRM**2)
         END IF
C
         IF (NUMINF.NE.0) THEN
C
C           Compute the transformed gradient of either the sum of
C           infeasibilities or the objective.  Initialize
C           SINGLR and UNITGZ.
C
            CALL E04NCP(PRBTYP,LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,
     *                  LDA,LDZY,LDR,NRANK,NZ,NRZ,ISTATE,KX,BIGBND,
     *                  TOLRNK,NUMINF,SUMINF,BL,BU,A,W(LRES),FEATOL,
     *                  W(LGQ),W(LCQ),R,X,W(LWTINF),W(LZY),W(LWRK))
            IF (NUMINF.EQ.0 .AND. PRBTYP.NE.'FP' .AND. NPHASE.EQ.1) THEN
               ITMAX = ITER + ITMAX2
               NPHASE = 2
            END IF
         END IF
      END IF
C
      GZNORM = ZERO
      IF (NZ.GT.0) GZNORM = DNRM2(NZ,W(LGQ),1)
C
      IF (NRZ.EQ.NZ) THEN
         GRZNRM = GZNORM
      ELSE
         GRZNRM = ZERO
         IF (NRZ.GT.0) GRZNRM = DNRM2(NRZ,W(LGQ),1)
      END IF
C
      GFNORM = GZNORM
      IF (NFREE.GT.0 .AND. NACTIV.GT.0) GFNORM = DNRM2(NFREE,W(LGQ),1)
C
C        ------------------------------------------------------------
C        Print the details of this iteration.
C        ------------------------------------------------------------
C        Define small quantities that reflect the magnitude of  x,
C        R  and  the matrix of constraints in the working set.
C        Use the largest and smallest diagonals of  R  to estimate
C        the condition number of  Rz1.
C        Note that NRZ <= NRANK + 1.
C
      IF (NRZ.EQ.0) THEN
         SINGLR = .FALSE.
      ELSE
         IF (NUMINF.GT.0 .OR. NRZ.GT.NRANK) THEN
            ABSRZZ = ZERO
            SINGLR = .TRUE.
         ELSE
            CALL F06FLF(NRZ,R,LDR+1,DRZMAX,DRZMIN)
            ABSRZZ = ABS(R(NRZ,NRZ))
            ROWNRM = DNRM2(N,R(1,1),LDR)
            SINGLR = ABSRZZ .LE. DRZMAX*TOLRNK .OR. ROWNRM .LE.
     *               TOLRNK .OR. ABS(R(1,1)) .LE. ROWNRM*TOLRNK
         END IF
C
      END IF
C
      CONDRZ = F06BLF(DRZMAX,DRZMIN,OVERFL)
C
      CONDT = ONE
      IF (NACTIV.GT.0) CONDT = F06BLF(DTMAX,DTMIN,OVERFL)
C
      IF (PRNT) THEN
         CALL E04NCJ(PRBTYP,ISDEL,ITER,JADD,JDEL,MSGLVL,NACTIV,NFREE,N,
     *               NCLIN,NRANK,LDR,LDT,NZ,NRZ,ISTATE,ALFA,CONDRZ,
     *               CONDT,GRZNRM,NUMINF,SUMINF,CTX,SSQ,AX,R,W(LT),X,
     *               W(LWRK))
C
         JDEL = 0
         JADD = 0
         ALFA = ZERO
      END IF
C
      IF (NUMINF.GT.0) THEN
         DINKY = ZERO
      ELSE
         OBJSIZ = ONE + ABS(SSQ+CTX)
         WSSIZE = ZERO
         IF (NACTIV.GT.0) WSSIZE = DTMAX
         DINKY = EPSPT8*MAX(WSSIZE,OBJSIZ,GFNORM)
         IF (UNCON) THEN
            UNITGZ = GRZNRM .LE. DINKY
         END IF
      END IF
C
C     If the projected gradient  Z'g  is small and Rz is of full
C     rank, X is a minimum on the working set.  An additional
C     refinement step is allowed to take care of an inaccurate
C     value of DINKY.
C
      STATPT = .NOT. SINGLR .AND. GRZNRM .LE. DINKY .OR. IREFN .GT.
     *         MREFN
C
      IF ( .NOT. STATPT) THEN
C        ---------------------------------------------------------
C        Compute a search direction.
C        ---------------------------------------------------------
         PRNT = .TRUE.
C
         ERROR = ITER .GE. ITMAX
         IF ( .NOT. ERROR) THEN
C
            IREFN = IREFN + 1
            ITER = ITER + 1
C
            CALL E04NCQ(LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,LDA,
     *                  LDZY,LDR,NRANK,NUMINF,NRZ,KX,CTP,PNORM,A,W(LAP),
     *                  W(LRES),W(LHZ),W(LPX),W(LGQ),W(LCQ),R,W(LZY),
     *                  W(LWRK))
C
C           ------------------------------------------------------
C           Find the constraint we bump into along P.
C           Update X and AX if the step ALFA is nonzero.
C           ------------------------------------------------------
C           ALFHIT is initialized to BIGALF.  If it remains
C           that way after the call to E04UCG, it will be
C           regarded as infinite.
C
            BIGALF = F06BLF(BIGDX,PNORM,OVERFL)
C
            CALL E04UCG(FIRSTV,HITLOW,ISTATE,INFORM,JADD,N,NCTOTL,
     *                  NUMINF,ALFHIT,PALFA,ATPHIT,BIGALF,BIGBND,PNORM,
     *                  W(LANORM),W(LAP),AX,BL,BU,FEATOL,W(LPX),X)
C
C           If  Rz1  is nonsingular,  ALFA = 1.0  will be the
C           step to the least-squares minimizer on the
C           current subspace. If the unit step does not violate
C           the nearest constraint by more than FEATOL,  the
C           constraint is not added to the working set.
C
            HITCON = SINGLR .OR. PALFA .LE. ONE
            UNCON = .NOT. HITCON
C
            IF (HITCON) THEN
               ALFA = ALFHIT
            ELSE
               JADD = 0
               ALFA = ONE
            END IF
C
C           Check for an unbounded solution or negligible step.
C
            UNBNDD = ALFA .GE. BIGALF
            STALL = ABS(ALFA*PNORM) .LE. EPSPT9*XNORM
            IF (STALL) THEN
               NSTALL = NSTALL + 1
               CYCLIN = NSTALL .GT. MSTALL
            ELSE
               NSTALL = 0
            END IF
C
            ERROR = UNBNDD .OR. CYCLIN
            IF ( .NOT. ERROR) THEN
C              ---------------------------------------------------
C              Set X = X + ALFA*P.  Update AX, GQ, RES and CTX.
C              ---------------------------------------------------
               IF (ALFA.NE.ZERO) CALL E04NCL(HITCON,HITLOW,LINOBJ,
     *                                UNITGZ,NCLIN,NRANK,NRZ,N,LDR,JADD,
     *                                NUMINF,ALFA,CTP,CTX,XNORM,W(LAP),
     *                                AX,BL,BU,W(LGQ),W(LHZ),W(LPX),
     *                                W(LRES),R,X,W(LWRK))
C
               IF (HITCON) THEN
C                 ------------------------------------------------
C                 Add a constraint to the working set.
C                 Update the TQ factors of the working set.
C                 Use P as temporary work space.
C                 ------------------------------------------------
C                 Update  ISTATE.
C
                  IF (BL(JADD).EQ.BU(JADD)) THEN
                     ISTATE(JADD) = 3
                  ELSE IF (HITLOW) THEN
                     ISTATE(JADD) = 1
                  ELSE
                     ISTATE(JADD) = 2
                  END IF
                  IADD = JADD - N
                  IF (JADD.LE.N) THEN
C
                     DO 40 IFIX = 1, NFREE
                        IF (KX(IFIX).EQ.JADD) GO TO 60
   40                CONTINUE
                  END IF
   60             CONTINUE
C
                  CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,
     *                        NFREE,NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,
     *                        KX,CONDMX,A,R,W(LT),W(LRES),W(LGQ),W(LZY),
     *                        W(LWRK),W(LRLAM),W(LPX),MSGLVL)
C
                  NRZ = NRZ - 1
                  NZ = NZ - 1
C
                  IF (JADD.LE.N) THEN
C
C                    A simple bound has been added.
C
                     NFREE = NFREE - 1
                  ELSE
C
C                    A general constraint has been added.
C
                     NACTIV = NACTIV + 1
                     KACTIV(NACTIV) = IADD
                  END IF
C
                  IREFN = 0
               END IF
C
C              ---------------------------------------------------
C              Check the feasibility of constraints with non-
C              negative ISTATE values.  If some violations have
C              occurred.  Refine the current X and set INFORM so
C              that feasibility is checked in E04NCP.
C              ---------------------------------------------------
               CALL E04NCR(N,NCLIN,ISTATE,BIGBND,CNORM,ERR1,JMAX1,NVIOL,
     *                     AX,BL,BU,FEATOL,X,W(LWRK))
C
               IF (ERR1.GT.FEATOL(JMAX1)) THEN
                  CALL E04NCH(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,
     *                        NRANK,NZ,N,NCTOTL,LDZY,LDA,LDR,LDT,ISTATE,
     *                        KACTIV,KX,JMAX1,ERR2,CTX,XNORM,A,AX,BL,BU,
     *                        W(LCQ),W(LRES),W(LRES0),FEATOL,R,W(LT),X,
     *                        W(LZY),W(LPX),W(LWRK))
C
                  IF (ROWERR) THEN
                     IF (MSGLVL.GT.0) THEN
                        WRITE (REC,FMT=99998)
                        CALL X04BAF(IPRINT,REC(1))
                     END IF
                     NUMINF = 1
                     ERROR = .TRUE.
                  ELSE
                     NUMINF = -1
                     UNCON = .FALSE.
                     IREFN = 0
                  END IF
               END IF
               NEEDFG = ALFA .NE. ZERO
            END IF
         END IF
      END IF
C
C        UNTIL      STATPT  .OR.  ERROR
      IF ( .NOT. (STATPT .OR. ERROR)) GO TO 20
C
C     ===============================================================
C     Try and find the index JDEL of a constraint to drop from
C     the working set.
C     ===============================================================
      JDEL = 0
C
      IF (NUMINF.EQ.0 .AND. PRBTYP.EQ.'FP') THEN
         IF (N.GT.NZ) CALL F06FBF(N-NZ,(ZERO),W(LRLAM),1)
         JTINY = 0
         JSMLST = 0
         JBIGST = 0
      ELSE
C
         CALL E04NCK(PRBTYP,MSGLVL,N,NACTIV,NFREE,LDA,LDT,NUMINF,NZ,NRZ,
     *               ISTATE,KACTIV,KX,DINKY,JSMLST,KSMLST,JINF,JTINY,
     *               JBIGST,KBIGST,TRULAM,A,W(LANORM),W(LGQ),W(LRLAM),
     *               W(LT),W(LWTINF))
C
      END IF
C
      IF ( .NOT. ERROR) THEN
         IF (JSMLST.GT.0) THEN
C
C           E04NCK found a regular constraint with multiplier less
C           than (-DINKY).
C
            JDEL = JSMLST
            KDEL = KSMLST
            ISDEL = ISTATE(JDEL)
            ISTATE(JDEL) = 0
C
         ELSE IF (JSMLST.LT.0) THEN
C
            JDEL = JSMLST
C
         ELSE IF (NUMINF.GT.0 .AND. JBIGST.GT.0) THEN
C
C           No feasible point exists for the constraints but the
C           sum of the constraint violations may be reduced by
C           moving off constraints with multipliers greater than 1.
C
            JDEL = JBIGST
            KDEL = KBIGST
            ISDEL = ISTATE(JDEL)
            IF (TRULAM.LE.ZERO) IS = -1
            IF (TRULAM.GT.ZERO) IS = -2
            ISTATE(JDEL) = IS
            FIRSTV = .TRUE.
            NUMINF = NUMINF + 1
         END IF
C
         IF (JDEL.NE.0 .AND. SINGLR) THEN
C
C           Cannot delete a constraint when Rz is singular.
C           Probably a weak minimum.
C
            JDEL = 0
         ELSE IF (JDEL.NE.0) THEN
C
C           Constraint JDEL has been deleted.
C           Update the matrix factorizations.
C
            CALL E04NCT(UNITQ,N,NACTIV,NFREE,NRES,NGQ,NZ,NRZ,LDA,LDZY,
     *                  LDR,LDT,NRANK,JDEL,KDEL,KACTIV,KX,A,W(LRES),R,
     *                  W(LT),W(LGQ),W(LZY),W(LWRK),W(LPX))
         END IF
      END IF
C
      IREFN = 0
      CONVRG = JDEL .EQ. 0
C
      PRNT = .FALSE.
      UNCON = .FALSE.
      NEEDFG = .FALSE.
C
C     until       convrg  .or.  error
      IF ( .NOT. (CONVRG .OR. ERROR)) GO TO 20
C
C     .....................End of main loop........................
C
      WEAK = JTINY .GT. 0 .OR. SINGLR
C
      IF (ERROR) THEN
         IF (UNBNDD) THEN
            INFORM = 2
            IF (NUMINF.GT.0) INFORM = 3
         ELSE IF (ITER.GE.ITMAX) THEN
            INFORM = 4
         ELSE IF (CYCLIN) THEN
            INFORM = 5
         END IF
      ELSE IF (CONVRG) THEN
         INFORM = 0
         IF (NUMINF.GT.0) THEN
            INFORM = 3
         ELSE IF (PRBTYP.NE.'FP' .AND. WEAK) THEN
            INFORM = 1
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Set   CLAMDA.  Print the full solution.
C     ------------------------------------------------------------------
      IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99999) PRBTYP, ITER
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      CALL E04UPP(NFREE,LDA,N,NCLIN,NCTOTL,NACTIV,ISTATE,KACTIV,KX,A,BL,
     *            BU,X,CLAMDA,FEATOL,W(LWRK),W(LRLAM),X)
      CALL E04MFK(MSGLVL,N,NCLIN,NCTOTL,BIGBND,NAMED,NAMES,ISTATE,BL,BU,
     *            CLAMDA,FEATOL,W(LWRK))
C
      RETURN
C
C
C     End of  E04NCZ. (LSCORE)
C
99999 FORMAT (/' Exit from ',A2,' problem after ',I5,' iterations.')
99998 FORMAT (' XXX  Warning.  Cannot satisfy the constraints to the a',
     *       'ccuracy requested.')
      END
