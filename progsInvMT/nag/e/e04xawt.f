      SUBROUTINE E04XAW(INFORM,LVLDER,MSGLVL,NCSET,N,NCNLN,LDCJ,LDCJU,
     *                  BIGBND,EPSRF,OKTOL,FDCHK,XNORM,CONFUN,NEEDC,BL,
     *                  BU,C,C1,CJAC,CJACU,CJDX,DX,ERR,X,Y,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-999 (JUN 1993).
C     MARK 17 REVISED. IER-1626 (JUN 1995).
C
C     ******************************************************************
C     E04XAW  checks if the gradients of the constraints have been coded
C     correctly.
C
C     On input,  the values of the constraints at the point X are stored
C     in C.  Their corresponding gradients are stored in CJACU.  If any
C     Jacobian element has not been specified,  it will have a dummy
C     value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves
C     satisfactory and no further information is desired, E04XAW is
C     terminated. Otherwise, E04XAZ is called to give optimal
C     step-sizes and a central-difference approximation to each
C     element of the Jacobian for which a test is deemed necessary,
C     either by the program or the user.
C
C     LVRFYC has the following meaning...
C
C     -1        do not perform any check.
C     0        do the cheap test only.
C     2 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     Level-2 matrix routines added 18-May-1988.
C     This version of E04XAW dated 14-Sep-1992.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
      DOUBLE PRECISION  ZERO, HALF, POINT9
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,POINT9=0.9D+0)
      DOUBLE PRECISION  ONE, TWO, TEN
      PARAMETER         (ONE=1.0D+0,TWO=2.0D+0,TEN=1.0D+1)
      CHARACTER*4       LBAD, LGOOD
      PARAMETER         (LBAD='BAD?',LGOOD='  OK')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, EPSRF, FDCHK, OKTOL, XNORM
      INTEGER           INFORM, LDCJ, LDCJU, LVLDER, MSGLVL, N, NCNLN,
     *                  NCSET
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), CJDX(*), DX(N), ERR(*), USER(*),
     *                  X(N), Y(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, LVRFYC, NOUT
C     .. Arrays in Common ..
      INTEGER           JVERFY(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CDEST, CIJ, CJDIFF, CJSIZE,
     *                  COLMAX, DXJ, DXMULT, EMAX, EPSACI, ERRBND, F1,
     *                  F2, FDEST, H, HOPT, HPHI, SDEST, SIGNH, STEPBL,
     *                  STEPBU, XJ
      INTEGER           I, IMAX, INFO, IROW, ITER, ITMAX, J, J3, J4,
     *                  JCOL, MODE, NCHECK, NCOLJ, NGOOD, NSTATE, NWRONG
      LOGICAL           CONST, DEBUG, DONE, FIRST, HEADNG, NEEDED, OK
      CHARACTER*4       KEY
C     .. Local Arrays ..
      CHARACTER*18      RESULT(0:4)
      CHARACTER*120     REC(4)
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, E04XAZ, F06DBF, F06FBF,
     *                  F06QFF, F06QHF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
C     .. Data statements ..
      DATA              RESULT/'                 ', 'Constant?      ',
     *                  'Linear or odd?   ', 'Too nonlinear?',
     *                  'Small derivative?'/
C     .. Executable Statements ..
C
      INFORM = 0
      NEEDED = NCNLN .GT. 0 .AND. LVRFYC .EQ. 0 .OR. LVRFYC .EQ. 2 .OR.
     *         LVRFYC .EQ. 3
      IF ( .NOT. NEEDED) RETURN
C
      IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,4,REC)
      END IF
      DEBUG = .FALSE.
      NSTATE = 0
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Perform the cheap test.
C     ==================================================================
      H = (ONE+XNORM)*FDCHK
C
      IF (N.LE.100) THEN
         DXMULT = 0.9D0
      ELSE IF (N.LE.250) THEN
         DXMULT = 0.99D0
      ELSE
         DXMULT = 0.999D0
      END IF
C
      DXJ = ONE/N
      DO 20 J = 1, N
         DX(J) = DXJ
         DXJ = -DXJ*DXMULT
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Do not perturb  X(J)  if the  J-th  column contains any
C     unknown elements.  Compute the directional derivative for each
C     constraint gradient.
C     ------------------------------------------------------------------
      NCHECK = 0
      DO 60 J = 1, N
         DO 40 I = 1, NCNLN
            IF (CJAC(I,J).EQ.RDUMMY) THEN
               DX(J) = ZERO
               GO TO 60
            END IF
   40    CONTINUE
         NCHECK = NCHECK + 1
C
         XJ = X(J)
         STEPBL = -ONE
         STEPBU = ONE
         IF (BL(J).GT.BIGLOW) STEPBL = MAX(STEPBL,BL(J)-XJ)
         IF (BU(J).LT.BIGUPP .AND. BU(J).GT.BL(J)) STEPBU = MIN(STEPBU,
     *       BU(J)-XJ)
C
         IF (HALF*(STEPBL+STEPBU).LT.ZERO) THEN
            DX(J) = DX(J)*STEPBL
         ELSE
            DX(J) = DX(J)*STEPBU
         END IF
   60 CONTINUE
C
      IF (NCHECK.EQ.0) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99995)
            CALL X04BAY(IPRINT,2,REC)
         END IF
      ELSE
C
C        Compute  (Jacobian)*DX.
C
         CALL DGEMV('Normal',NCNLN,N,ONE,CJACU,LDCJU,DX,1,ZERO,CJDX,1)
C
C        ---------------------------------------------------------------
C        Make forward-difference approximation along DX.
C        ---------------------------------------------------------------
         CALL DCOPY(N,X,1,Y,1)
         CALL DAXPY(N,H,DX,1,Y,1)
C
         CALL F06DBF(NCNLN,(1),NEEDC,1)
C
         MODE = 0
         CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,Y,C1,CJACU,NSTATE,IUSER,
     *               USER)
         IF (MODE.LT.0) GO TO 160
C
C        Set  ERR = (C1 - C)/H  - Jacobian*DX.  This should be small.
C
         DO 80 I = 1, NCNLN
            ERR(I) = (C1(I)-C(I))/H - CJDX(I)
   80    CONTINUE
         IMAX = IDAMAX(NCNLN,ERR,1)
         EMAX = ABS(ERR(IMAX))/(ABS(CJDX(IMAX))+ONE)
C
         IF (MSGLVL.GT.0) THEN
            IF (EMAX.LE.OKTOL) THEN
               WRITE (REC,FMT=99998)
               CALL X04BAY(IPRINT,2,REC)
            ELSE
               WRITE (REC,FMT=99997)
               CALL X04BAY(IPRINT,2,REC)
            END IF
            WRITE (REC,FMT=99996) EMAX, IMAX
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (EMAX.GE.POINT9) INFORM = 1
      END IF
C
C     ==================================================================
C     Element-wise check.
C     ==================================================================
      IF (LVRFYC.GE.2) THEN
         IF (LVLDER.EQ.3) THEN
C
C           Recompute the Jacobian to find the non-constant elements.
C
            CALL F06QHF('General',NCNLN,N,RDUMMY,RDUMMY,CJACU,LDCJU)
C
            CALL F06DBF(NCNLN,(1),NEEDC,1)
            NSTATE = 0
            MODE = 2
C
            CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,NSTATE,
     *                  IUSER,USER)
            IF (MODE.LT.0) GO TO 160
C
         END IF
C
         CALL F06DBF(NCNLN,(0),NEEDC,1)
C
         ITMAX = 3
         NCHECK = 0
         NWRONG = 0
         NGOOD = 0
         COLMAX = -ONE
         JCOL = 0
         IROW = 0
         MODE = 0
         J3 = JVERFY(3)
         J4 = JVERFY(4)
C
C        ---------------------------------------------------------------
C        Loop over each column.
C        ---------------------------------------------------------------
         DO 140 J = J3, J4
C
            CALL F06FBF(NCNLN,ZERO,ERR,1)
            NCOLJ = 0
            HEADNG = .TRUE.
            XJ = X(J)
C
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
            SIGNH = ONE
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) SIGNH = -ONE
C
            DO 120 I = 1, NCNLN
               EPSACI = EPSRF*(ONE+ABS(C(I)))
C
               IF (CJACU(I,J).NE.RDUMMY) THEN
C                 ------------------------------------------------------
C                 Check this Jacobian element.
C                 ------------------------------------------------------
                  NCHECK = NCHECK + 1
                  NCOLJ = NCOLJ + 1
                  NEEDC(I) = 1
C
                  CIJ = CJAC(I,J)
                  CJSIZE = ABS(CIJ)
C                 ------------------------------------------------------
C                 Find a finite-difference interval by iteration.
C                 ------------------------------------------------------
                  ITER = 0
                  HOPT = TWO*(ONE+ABS(XJ))*SQRT(EPSRF)
                  H = TEN*HOPT*SIGNH
                  CDEST = ZERO
                  SDEST = ZERO
                  FIRST = .TRUE.
C
C                 +                REPEAT
  100             X(J) = XJ + H
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 160
                  F1 = C1(I)
C
                  X(J) = XJ + H + H
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 160
                  F2 = C1(I)
C
                  CALL E04XAZ(DEBUG,DONE,FIRST,EPSACI,EPSRF,C(I),INFO,
     *                        ITER,ITMAX,CDEST,FDEST,SDEST,ERRBND,F1,F2,
     *                        H,HOPT,HPHI)
C
C                 +                UNTIL     DONE
                  IF ( .NOT. DONE) GO TO 100
C
C                 ------------------------------------------------------
C                 Exit for this element.
C                 ------------------------------------------------------
                  CJDIFF = CDEST
                  ERR(I) = ABS(CJDIFF-CIJ)/(CJSIZE+ONE)
C
                  OK = ERR(I) .LE. OKTOL
                  IF (OK) THEN
                     KEY = LGOOD
                     NGOOD = NGOOD + 1
                  ELSE
                     KEY = LBAD
                     NWRONG = NWRONG + 1
                  END IF
C
                  IF (MSGLVL.GT.0) THEN
                     CONST = OK .AND. INFO .EQ. 1 .AND. ABS(CIJ)
     *                       .LT. EPSPT8
                     IF ( .NOT. CONST) THEN
                        IF (HEADNG) THEN
                           WRITE (REC,FMT=99994)
                           CALL X04BAY(IPRINT,4,REC)
                           IF (OK) THEN
                              WRITE (REC,FMT=99993) J, XJ, HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER
                           ELSE
                              WRITE (REC,FMT=99992) J, XJ, HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER, RESULT(INFO)
                           END IF
                           CALL X04BAF(IPRINT,REC(1))
                           HEADNG = .FALSE.
                        ELSE
                           IF (OK) THEN
                              WRITE (REC,FMT=99991) HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER
                           ELSE
                              WRITE (REC,FMT=99990) HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER, RESULT(INFO)
                           END IF
                           CALL X04BAF(IPRINT,REC(1))
                        END IF
                     END IF
                  END IF
                  NEEDC(I) = 0
               END IF
  120       CONTINUE
C
C           ------------------------------------------------------------
C           Finished with this column.
C           ------------------------------------------------------------
            IF (NCOLJ.GT.0) THEN
               IMAX = IDAMAX(NCNLN,ERR,1)
               EMAX = ABS(ERR(IMAX))
C
               IF (EMAX.GE.COLMAX) THEN
                  IROW = IMAX
                  JCOL = J
                  COLMAX = EMAX
               END IF
            END IF
            X(J) = XJ
C
  140    CONTINUE
C
         INFORM = 0
         IF (COLMAX.GE.POINT9) INFORM = 1
C
         IF (MSGLVL.GT.0) THEN
            IF (NCHECK.EQ.0) THEN
               WRITE (REC,FMT=99986) NCSET
               CALL X04BAF(IPRINT,REC(1))
            ELSE
               IF (NWRONG.EQ.0) THEN
                  WRITE (REC,FMT=99989) NGOOD, NCHECK, J3, J4
               ELSE
                  WRITE (REC,FMT=99988) NWRONG, NCHECK, J3, J4
               END IF
               CALL X04BAY(IPRINT,3,REC)
               WRITE (REC,FMT=99987) COLMAX, IROW, JCOL
               CALL X04BAY(IPRINT,3,REC)
            END IF
         END IF
C
      END IF
C
C     Copy  ( constants + gradients + dummy values )  back into CJACU.
C
      CALL F06QFF('General',NCNLN,N,CJAC,LDCJ,CJACU,LDCJU)
C
      RETURN
C
  160 INFORM = MODE
      RETURN
C
C
C     End of  E04XAW. (CHCJAC)
C
99999 FORMAT (//' Verification of the constraint gradients.',/' ------',
     *       '-----------------------------------')
99998 FORMAT (/' The constraint Jacobian seems to be ok.')
99997 FORMAT (/' XXX  The constraint Jacobian seems to be incorrect.')
99996 FORMAT (/' The largest relative error was',1P,D12.2,'  in constr',
     *       'aint',I5,/)
99995 FORMAT (/' Every column contains a constant or missing element.')
99994 FORMAT (//' Column    X(J)     DX(J)    Row    Jacobian Value   ',
     *       '   Difference Approxn  Itns',/)
99993 FORMAT (I7,1P,2D10.2,I5,1P,2D18.8,2X,A4,I6)
99992 FORMAT (I7,1P,2D10.2,I5,1P,2D18.8,2X,A4,I6,2X,A18)
99991 FORMAT (17X,1P,D10.2,I5,1P,2D18.8,2X,A4,I6)
99990 FORMAT (17X,1P,D10.2,I5,1P,2D18.8,2X,A4,I6,2X,A18)
99989 FORMAT (/I7,'  constraint Jacobian elements out of the',I6,/9X,
     *       'set in cols',I6,'  through',I6,'  seem to be ok.')
99988 FORMAT (/' XXX  There seem to be',I6,'  incorrect Jacobian eleme',
     *       'nts out of the',I6,/8X,'set in cols',I6,'  through',I6)
99987 FORMAT (/' The largest relative error was',1P,D12.2,'  in row',I5,
     *       ',  column',I5,/)
99986 FORMAT (' All',I6,'   assigned Jacobian elements are constant.')
      END
