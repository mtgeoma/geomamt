      SUBROUTINE E04UPR(INFORM,MSGLVL,LVLDER,M,N,NCNLN,LDCJ,LDCJU,LDFJ,
     *                  LDFJU,BIGBND,EPSRF,FDNORM,CONFUN,OBJFUN,NEEDC,
     *                  BL,BU,C,C1,C2,CJAC,CJACU,F,F1,F2,FJAC,FJACU,
     *                  HFORWD,HCNTRL,X,Y,IUSER,USER)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1098 (JUL 1993).
C     MARK 17 REVISED. IER-1619 (JUN 1995).
C
C     ******************************************************************
C     E04UPR  computes difference intervals for missing Jacobian
C     elements of f(x) and c(x).  Intervals are computed using a
C     procedure that requires about two function evaluations if the
C     function is well scaled.  Central-difference gradients are
C     obtained as a by-product of the algorithm.
C
C     On entry...
C        F and C contain the problem functions at the point X.
C        An element of CJAC or FJAC not equal to RDUMMY signifies a know
C        gradient value.  Such values are not estimated by differencing.
C        CJACU and FJACU have dummy elements in the same positions as
C        CJAC and FJAC.
C
C     On exit...
C        CJAC and FJAC contain central-difference derivative estimates.
C        Elements of CJACU and FJACU are unaltered except for those
C        corresponding to constant derivatives, which are given the same
C        values as CJAC or FJAC.
C
C     Systems Optimization Laboratory,
C     Department of Operations Research,
C     Stanford University, Stanford, California 94305
C     Original version written 10-May-1988.
C     This version of E04UPR dated 14-Sep-1992.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
      DOUBLE PRECISION  FACTOR
      PARAMETER         (FACTOR=0.97D+0)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      DOUBLE PRECISION  TWO, FOUR, TEN
      PARAMETER         (TWO=2.0D+0,FOUR=4.0D+0,TEN=1.0D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, EPSRF, FDNORM
      INTEGER           INFORM, LDCJ, LDCJU, LDFJ, LDFJU, LVLDER, M,
     *                  MSGLVL, N, NCNLN
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), C2(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), F(M), F1(M), F2(M),
     *                  FJAC(LDFJ,*), FJACU(LDFJU,*), HCNTRL(*),
     *                  HFORWD(*), USER(*), X(N), Y(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LFDSET, LINES1, LINES2, LVLDIF,
     *                  NCDIFF, NFDIFF, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CDEST, CJDIFF, D, DX, EPSA,
     *                  ERRBND, ERRMAX, ERRMIN, FBACK, FDEST, FFORW,
     *                  FJDIFF, FX, H, HCD, HFD, HMAX, HMIN, HOPT, HPHI,
     *                  SDEST, SIGNH, STEPBL, STEPBU, SUMEPS, SUMSD,
     *                  TEST, XJ, YJ
      INTEGER           I, INFO, IROW, IROW1, IROW2, ITER, ITMAX, J,
     *                  MODE, NCCNST, NCOLJ, NFCNST, NSTATE
      LOGICAL           DEBUG, DONE, FIRST, HEADNG, NEEDED
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          E04XAZ, F06DBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Executable Statements ..
C
      INFORM = 0
      NEEDED = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2 .OR. LVLDER .EQ. 1 .AND.
     *         NCNLN .GT. 0
      IF ( .NOT. NEEDED) RETURN
C
      DEBUG = .FALSE.
      IF (LFDSET.EQ.0) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99999)
            CALL X04BAY(IPRINT,4,REC)
         END IF
C
         NSTATE = 0
         ITMAX = 3
         MODE = 0
C
         NCCNST = 0
         NFCNST = 0
         HEADNG = .TRUE.
C
         FDNORM = ZERO
C
C        ===============================================================
C        For each column of the matrix
C              ( CJAC  )
C              ( FJAC  ),
C        rows IROW1 thru IROW2 are searched for missing elements.
C        ===============================================================
         IROW1 = 1
         IROW2 = NCNLN + M
         IF (LVLDER.EQ.1) IROW2 = NCNLN
         IF (LVLDER.EQ.2) IROW1 = NCNLN + 1
C
         BIGLOW = -BIGBND
         BIGUPP = BIGBND
C
         IF (NCNLN.GT.0) CALL F06DBF(NCNLN,(0),NEEDC,1)
C
         DO 60 J = 1, N
            XJ = X(J)
            NCOLJ = 0
            SUMSD = ZERO
            SUMEPS = ZERO
            HFD = ZERO
            HCD = ZERO
            HMAX = ZERO
            HMIN = ONE/EPSPT3
            ERRMAX = ZERO
            ERRMIN = ZERO
C
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
            SIGNH = ONE
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) SIGNH = -ONE
C
            DO 40 IROW = IROW1, IROW2
C
               IF (IROW.LE.NCNLN) THEN
                  I = IROW
                  TEST = CJACU(I,J)
               ELSE
                  I = IROW - NCNLN
                  TEST = FJACU(I,J)
               END IF
C
               IF (TEST.EQ.RDUMMY) THEN
C                 ======================================================
C                 Get the difference interval for this element.
C                 ======================================================
                  NCOLJ = NCOLJ + 1
C
                  IF (IROW.LE.NCNLN) THEN
                     NEEDC(I) = 1
                     FX = C(I)
                  ELSE
                     FX = F(I)
                  END IF
C
                  EPSA = EPSRF*(ONE+ABS(FX))
C
C                 ------------------------------------------------------
C                 Find a finite-difference interval by iteration.
C                 ------------------------------------------------------
                  ITER = 0
                  HOPT = TWO*(ONE+ABS(XJ))*SQRT(EPSRF)
                  H = SIGNH*TEN*HOPT
                  CDEST = ZERO
                  SDEST = ZERO
                  FIRST = .TRUE.
C
C                 +                REPEAT
   20             X(J) = XJ + H
                  IF (IROW.LE.NCNLN) THEN
                     CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                           NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 240
                     FFORW = C1(I)
                  ELSE
                     CALL OBJFUN(MODE,M,N,LDFJU,X,F1,FJACU,NSTATE,IUSER,
     *                           USER)
                     IF (MODE.LT.0) GO TO 240
                     FFORW = F1(I)
                  END IF
C
                  X(J) = XJ + H + H
                  IF (IROW.LE.NCNLN) THEN
                     CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                           NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 240
                     FBACK = C1(I)
                  ELSE
                     CALL OBJFUN(MODE,M,N,LDFJU,X,F1,FJACU,NSTATE,IUSER,
     *                           USER)
                     IF (MODE.LT.0) GO TO 240
                     FBACK = F1(I)
                  END IF
C
                  CALL E04XAZ(DEBUG,DONE,FIRST,EPSA,EPSRF,FX,INFO,ITER,
     *                        ITMAX,CDEST,FDEST,SDEST,ERRBND,FFORW,
     *                        FBACK,H,HOPT,HPHI)
C
C                 +                UNTIL     DONE
                  IF ( .NOT. DONE) GO TO 20
C
                  IF (IROW.LE.NCNLN) THEN
                     CJAC(I,J) = CDEST
                     IF (INFO.EQ.1 .OR. INFO.EQ.2) THEN
                        NCCNST = NCCNST + 1
                        NCDIFF = NCDIFF - 1
                        CJACU(I,J) = -RDUMMY
                     END IF
                  ELSE
                     FJAC(I,J) = CDEST
                     IF (INFO.EQ.1 .OR. INFO.EQ.2) THEN
                        NFCNST = NFCNST + 1
                        NFDIFF = NFDIFF - 1
                        FJACU(I,J) = -RDUMMY
                     END IF
                  END IF
C
                  SUMSD = SUMSD + ABS(SDEST)
                  SUMEPS = SUMEPS + EPSA
                  IF (HOPT.GT.HMAX) THEN
                     HMAX = HOPT
                     ERRMAX = ERRBND
                  END IF
                  IF (HOPT.LT.HMIN) THEN
                     HMIN = HOPT
                     ERRMIN = ERRBND
                  END IF
C
                  IF (INFO.EQ.0) HCD = MAX(HCD,HPHI)
               END IF
   40       CONTINUE
C
            IF (NCOLJ.GT.0) THEN
               IF (HMIN.GT.HMAX) THEN
                  HMIN = HMAX
                  ERRMIN = ERRMAX
               END IF
C
               IF (FOUR*SUMEPS.LT.HMIN*HMIN*SUMSD) THEN
                  HFD = HMIN
                  ERRMAX = ERRMIN
               ELSE IF (FOUR*SUMEPS.GT.HMAX*HMAX*SUMSD) THEN
                  HFD = HMAX
               ELSE
                  HFD = TWO*SQRT(SUMEPS/SUMSD)
                  ERRMAX = TWO*SQRT(SUMEPS*SUMSD)
               END IF
C
               IF (HCD.EQ.ZERO) HCD = TEN*HFD
C
               IF (MSGLVL.GT.0) THEN
                  IF (HEADNG) THEN
                     WRITE (REC,FMT=99998)
                     CALL X04BAY(IPRINT,4,REC)
                  END IF
                  WRITE (REC,FMT=99997) J, XJ, HFD, HCD, ERRMAX
                  CALL X04BAF(IPRINT,REC(1))
                  HEADNG = .FALSE.
               END IF
               FDNORM = MAX(FDNORM,HFD)
               HFORWD(J) = HFD/(ONE+ABS(XJ))
               HCNTRL(J) = HCD/(ONE+ABS(XJ))
            END IF
            X(J) = XJ
   60    CONTINUE
C
         IF (NCCNST+NFCNST.GT.0) THEN
C
C           Check that the constants have been set properly by
C           evaluating the gradients at a strange (but feasible) point.
C
            D = ONE/N
C
            DO 80 J = 1, N
               XJ = X(J)
               STEPBL = -ONE
               STEPBU = ONE
               IF (BL(J).GT.BIGLOW) STEPBL = MAX(STEPBL,BL(J)-XJ)
               IF (BU(J).LT.BIGUPP .AND. BU(J).GT.BL(J))
     *             STEPBU = MIN(STEPBU,BU(J)-XJ)
C
               IF (HALF*(STEPBL+STEPBU).LT.ZERO) THEN
                  Y(J) = XJ + D*STEPBL
               ELSE
                  Y(J) = XJ + D*STEPBU
               END IF
C
               D = FACTOR*D
   80       CONTINUE
C
            IF (NCNLN.GT.0) THEN
               CALL F06DBF(NCNLN,(1),NEEDC,1)
               CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,Y,C2,CJACU,NSTATE,
     *                     IUSER,USER)
               IF (MODE.LT.0) GO TO 240
            END IF
C
            IF (M.GT.0) THEN
               CALL OBJFUN(MODE,M,N,LDFJU,Y,F2,FJACU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 240
            END IF
C
C           ------------------------------------------------------------
C           Loop over each of the elements of  x.
C           ------------------------------------------------------------
            DO 180 J = 1, N
               YJ = Y(J)
               DX = HALF*(X(J)-YJ)
               Y(J) = YJ + DX
C
               IF (NCNLN.GT.0) THEN
                  NCOLJ = 0
                  DO 100 I = 1, NCNLN
                     IF (CJACU(I,J).EQ.-RDUMMY) THEN
                        NEEDC(I) = 1
                        NCOLJ = NCOLJ + 1
                     ELSE
                        NEEDC(I) = 0
                     END IF
  100             CONTINUE
C
                  IF (NCOLJ.GT.0) THEN
                     CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,Y,C1,CJACU,
     *                           NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 240
C
                     DO 120 I = 1, NCNLN
                        IF (NEEDC(I).EQ.1) THEN
                           CJDIFF = (C1(I)-C2(I))/DX
                           IF (CJDIFF.EQ.CJAC(I,J)) THEN
                              CJACU(I,J) = CJDIFF
                           ELSE
                              CJACU(I,J) = RDUMMY
                              NCCNST = NCCNST - 1
                              NCDIFF = NCDIFF + 1
                           END IF
                        END IF
  120                CONTINUE
                  END IF
               END IF
C
C              Now check the objective Jacobian.
C
               NCOLJ = 0
               DO 140 I = 1, M
                  IF (FJACU(I,J).EQ.-RDUMMY) THEN
                     NCOLJ = NCOLJ + 1
                  END IF
  140          CONTINUE
C
               IF (NCOLJ.GT.0) THEN
                  CALL OBJFUN(MODE,M,N,LDFJU,Y,F1,FJACU,NSTATE,IUSER,
     *                        USER)
                  IF (MODE.LT.0) GO TO 240
C
                  DO 160 I = 1, M
                     IF (FJACU(I,J).EQ.-RDUMMY) THEN
                        FJDIFF = (F1(I)-F2(I))/DX
                        IF (FJDIFF.EQ.FJAC(I,J)) THEN
                           FJACU(I,J) = FJDIFF
                        ELSE
                           FJACU(I,J) = RDUMMY
                           NFCNST = NFCNST - 1
                           NFDIFF = NFDIFF + 1
                        END IF
                     END IF
  160             CONTINUE
               END IF
C
               Y(J) = YJ
  180       CONTINUE
C
            IF (MSGLVL.GT.0) THEN
               IF (LVLDER.LT.2 .AND. NCCNST.GT.0) THEN
                  WRITE (REC,FMT=99996) NCCNST
                  CALL X04BAY(IPRINT,2,REC)
               END IF
               IF (LVLDER.NE.1 .AND. NFCNST.GT.0) THEN
                  WRITE (REC,FMT=99995) NFCNST
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            END IF
C
            IF (NCDIFF.EQ.0 .AND. LVLDER.LT.2 .AND. NCNLN.GT.0) THEN
               IF (LVLDER.EQ.0) LVLDER = 2
               IF (LVLDER.EQ.1) LVLDER = 3
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99994) LVLDER
                  CALL X04BAY(IPRINT,4,REC)
               END IF
            END IF
C
            IF (NFDIFF.EQ.0 .AND. LVLDER.NE.1) THEN
               IF (LVLDER.EQ.0) LVLDER = 1
               IF (LVLDER.EQ.2) LVLDER = 3
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99993) LVLDER
                  CALL X04BAY(IPRINT,4,REC)
               END IF
            END IF
         END IF
      ELSE IF (LFDSET.EQ.2) THEN
C
C        The user has supplied HFORWD and HCNTRL.
C        Check for wild values.
C
         DO 200 J = 1, N
            IF (HFORWD(J).LE.ZERO) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99992) J, HFORWD(J), EPSPT5
                  CALL X04BAF(IPRINT,REC(1))
               END IF
               HFORWD(J) = EPSPT5
            END IF
  200    CONTINUE
         DO 220 J = 1, N
            IF (HCNTRL(J).LE.ZERO) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99991) J, HCNTRL(J), EPSPT3
                  CALL X04BAF(IPRINT,REC(1))
               END IF
               HCNTRL(J) = EPSPT3
            END IF
  220    CONTINUE
      END IF
C
      RETURN
C
  240 INFORM = MODE
      RETURN
C
C
C     End of  E04UPR. (CHFDLS)
C
99999 FORMAT (//' Computation of the finite-difference intervals',/' -',
     *       '---------------------------------------------')
99998 FORMAT (//'    J      X(J)   Forward DX(J)   Central DX(J)      ',
     *       'Error est.',/)
99997 FORMAT (I5,1P,D10.2,1P,D16.6,1P,2D16.6)
99996 FORMAT (/I5,'  constant constraint gradient elements assigned.')
99995 FORMAT (/I5,'  constant objective  gradient elements assigned.')
99994 FORMAT (//' All missing constraint Jacobian elements are constan',
     *       'ts.',/' Derivative level increased to ',I4)
99993 FORMAT (//' All missing objective  Jacobian elements are constan',
     *       'ts.',/' Derivative level increased to ',I4)
99992 FORMAT (' XXX  ',I4,'-th difference interval ',1P,D10.2,' replac',
     *       'ed by ',1P,D10.2)
99991 FORMAT (' XXX  ',I4,'-th central-difference interval ',1P,D10.2,
     *       ' replaced by ',1P,D10.2)
      END
