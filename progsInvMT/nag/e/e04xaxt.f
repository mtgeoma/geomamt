      SUBROUTINE E04XAX(INFORM,MSGLVL,N,BIGBND,EPSRF,OKTOL,FDCHK,OBJF,
     *                  XNORM,OBJFUN,BL,BU,GRAD,GRADU,DX,X,Y,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-882 (NOV 1990).
C     MARK 16A REVISED. IER-1000 (JUN 1993).
C     MARK 17 REVISED. IER-1627 (JUN 1995).
C
C     ******************************************************************
C     E04XAX  checks if the gradients of the objective function have
C     been coded correctly.
C
C     On input,  the value of the objective function at the point X is
C     stored in OBJF.  The corresponding gradient is stored in GRADU.
C     If any gradient element has not been specified,  it will have a
C     dummy value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods. If this proves
C     satisfactory and no further information is desired, E04XAX is
C     terminated. Otherwise, the routine E04XAZ is called to give
C     optimal step-sizes and a forward-difference approximation to
C     each element of the gradient for which a test is deemed
C     necessary, either by the program or the user.
C
C     Other inputs:
C
C        X         The n-dimensional point at which the
C                  gradient is to be verified.
C        EPSRF     The positive bound on the relative error
C                  associated with computing the function at
C                  the point x.
C        OKTOL     The desired relative accuracy which the
C                  elements of the gradient should satisfy.
C
C     LVRFYC has the following meaning...
C
C     -1        do not perform any check.
C     0        do the cheap test only.
C     1 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     This version of E04XAX  dated  19-Oct-1992.
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
      DOUBLE PRECISION  BIGBND, EPSRF, FDCHK, OBJF, OKTOL, XNORM
      INTEGER           INFORM, MSGLVL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), DX(N), GRAD(N), GRADU(N), USER(*),
     *                  X(N), Y(N)
      INTEGER           IUSER(*)
C     .. Subroutine Arguments ..
      EXTERNAL          OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, LVRFYC, NOUT
C     .. Arrays in Common ..
      INTEGER           JVERFY(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CDEST, DXJ, DXMULT, EMAX, EPSA,
     *                  ERRBND, ERROR, F1, F2, FDEST, GDIFF, GDX, GJ,
     *                  GSIZE, H, HOPT, HPHI, OBJF1, SDEST, STEPBL,
     *                  STEPBU, XJ
      INTEGER           INFO, ITER, ITMAX, J, J1, J2, JMAX, MODE,
     *                  NCHECK, NGOOD, NSTATE, NWRONG
      LOGICAL           CONST, DEBUG, DONE, FIRST, HEADNG, NEEDED, OK
      CHARACTER*4       KEY
C     .. Local Arrays ..
      CHARACTER*18      RESULT(0:4)
      CHARACTER*120     REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, E04XAZ, X04BAF, X04BAY
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
      NEEDED = LVRFYC .EQ. 0 .OR. LVRFYC .EQ. 1 .OR. LVRFYC .EQ. 3
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
C     Do not perturb X(J) if the  J-th  element is missing.
C     Compute the directional derivative.
C     ------------------------------------------------------------------
      NCHECK = 0
      DO 40 J = 1, N
         IF (GRAD(J).EQ.RDUMMY) THEN
            DX(J) = ZERO
         ELSE
            NCHECK = NCHECK + 1
C
            XJ = X(J)
            STEPBL = -ONE
            STEPBU = ONE
            IF (BL(J).GT.BIGLOW) STEPBL = MAX(STEPBL,BL(J)-XJ)
            IF (BU(J).LT.BIGUPP .AND. BU(J).GT.BL(J))
     *          STEPBU = MIN(STEPBU,BU(J)-XJ)
C
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) THEN
               DX(J) = DX(J)*STEPBL
            ELSE
               DX(J) = DX(J)*STEPBU
            END IF
         END IF
   40 CONTINUE
C
      IF (NCHECK.EQ.0) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99989)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         RETURN
      END IF
      GDX = DDOT(N,GRADU,1,DX,1)
C
C     ------------------------------------------------------------------
C     Make forward-difference approximation along  p.
C     ------------------------------------------------------------------
      CALL DCOPY(N,X,1,Y,1)
      CALL DAXPY(N,H,DX,1,Y,1)
C
      MODE = 0
      CALL OBJFUN(MODE,N,Y,OBJF1,GRADU,NSTATE,IUSER,USER)
      IF (MODE.LT.0) GO TO 100
C
      GDIFF = (OBJF1-OBJF)/H
      ERROR = ABS(GDIFF-GDX)/(ABS(GDX)+ONE)
C
      OK = ERROR .LE. OKTOL
C
      IF (MSGLVL.GT.0) THEN
         IF (OK) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAY(IPRINT,2,REC)
         ELSE
            WRITE (REC,FMT=99997)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         WRITE (REC,FMT=99996) GDX, GDIFF
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      IF (ERROR.GE.POINT9) INFORM = 1
C
C     ==================================================================
C     Element-wise check.
C     ==================================================================
      IF (LVRFYC.EQ.1 .OR. LVRFYC.EQ.3) THEN
         HEADNG = .TRUE.
         ITMAX = 3
         NCHECK = 0
         NWRONG = 0
         NGOOD = 0
         JMAX = 0
         EMAX = ZERO
         J1 = JVERFY(1)
         J2 = JVERFY(2)
C
C        ---------------------------------------------------------------
C        Loop over each of the elements of  x.
C        ---------------------------------------------------------------
         DO 80 J = J1, J2
C
            IF (GRAD(J).NE.RDUMMY) THEN
C              ---------------------------------------------------------
C              Check this gradient element.
C              ---------------------------------------------------------
               NCHECK = NCHECK + 1
               GJ = GRAD(J)
               GSIZE = ABS(GJ)
               XJ = X(J)
C              ---------------------------------------------------------
C              Find a finite-difference interval by iteration.
C              ---------------------------------------------------------
               ITER = 0
               EPSA = EPSRF*(ONE+ABS(OBJF))
               CDEST = ZERO
               SDEST = ZERO
               FIRST = .TRUE.
C
               STEPBL = BIGLOW
               STEPBU = BIGUPP
               IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
               IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
               HOPT = TWO*(ONE+ABS(XJ))*SQRT(EPSRF)
               H = TEN*HOPT
               IF (HALF*(STEPBL+STEPBU).LT.ZERO) H = -H
C
C              +             REPEAT
   60          X(J) = XJ + H
               CALL OBJFUN(MODE,N,X,F1,GRADU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 100
C
               X(J) = XJ + H + H
               CALL OBJFUN(MODE,N,X,F2,GRADU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 100
C
               CALL E04XAZ(DEBUG,DONE,FIRST,EPSA,EPSRF,OBJF,INFO,ITER,
     *                     ITMAX,CDEST,FDEST,SDEST,ERRBND,F1,F2,H,HOPT,
     *                     HPHI)
C
C              +             UNTIL     DONE
               IF ( .NOT. DONE) GO TO 60
C
C              ---------------------------------------------------------
C              Exit for this variable.
C              ---------------------------------------------------------
               GDIFF = CDEST
               X(J) = XJ
C
               ERROR = ABS(GDIFF-GJ)/(GSIZE+ONE)
               IF (ERROR.GE.EMAX) THEN
                  EMAX = ERROR
                  JMAX = J
               END IF
C
               OK = ERROR .LE. OKTOL
               IF (OK) THEN
                  KEY = LGOOD
                  NGOOD = NGOOD + 1
               ELSE
                  KEY = LBAD
                  NWRONG = NWRONG + 1
               END IF
C
               IF (MSGLVL.GT.0) THEN
C
C                 Zero elements are not printed.
C
                  CONST = OK .AND. INFO .EQ. 1 .AND. ABS(GJ) .LT. EPSPT8
                  IF ( .NOT. CONST) THEN
                     IF (HEADNG) THEN
                        WRITE (REC,FMT=99995)
                        CALL X04BAY(IPRINT,4,REC)
                        HEADNG = .FALSE.
                     END IF
                     IF (OK) THEN
                        WRITE (REC,FMT=99994) J, XJ, HOPT, GJ, GDIFF,
     *                    KEY, ITER
                     ELSE
                        WRITE (REC,FMT=99993) J, XJ, HOPT, GJ, GDIFF,
     *                    KEY, ITER, RESULT(INFO)
                     END IF
                     CALL X04BAF(IPRINT,REC(1))
                  END IF
               END IF
            END IF
   80    CONTINUE
C
C        ===============================================================
C        Done.
C        ===============================================================
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            IF (NWRONG.EQ.0) THEN
               WRITE (REC,FMT=99992) NGOOD, NCHECK, J1, J2
               CALL X04BAY(IPRINT,3,REC)
            ELSE
               WRITE (REC,FMT=99991) NWRONG, NCHECK, J1, J2
               CALL X04BAY(IPRINT,3,REC)
            END IF
            WRITE (REC,FMT=99990) EMAX, JMAX
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (ERROR.GE.POINT9) INFORM = 1
      END IF
C
      CALL DCOPY(N,GRAD,1,GRADU,1)
C
      RETURN
C
  100 INFORM = MODE
      RETURN
C
C
C     End of  E04XAX. (CHKGRD)
C
99999 FORMAT (//' Verification of the objective gradients.',/' -------',
     *       '---------------------------------')
99998 FORMAT (/' The objective gradients seem to be ok.')
99997 FORMAT (/' XXX  The objective gradients seem to be incorrect.')
99996 FORMAT (/' Directional derivative of the objective',1P,D18.8,
     *       /' Difference approximation               ',1P,D18.8)
99995 FORMAT (//4X,'J',4X,'X(J)',5X,'DX(J)',11X,'G(J)',11X,'Difference',
     *       ' approxn  Itns',/)
99994 FORMAT (I5,1P,2D10.2,1P,2D18.8,2X,A4,I6)
99993 FORMAT (I5,1P,2D10.2,1P,2D18.8,2X,A4,I6,2X,A18)
99992 FORMAT (/I7,'  Objective gradients out of the',I6,/9X,'set in co',
     *       'ls',I6,'  through',I6,'  seem to be ok.')
99991 FORMAT (/' XXX  There seem to be',I6,'  incorrect objective grad',
     *       'ients out of the',I6,/8X,'set in cols',I6,'  through',I6)
99990 FORMAT (/' The largest relative error was',1P,D12.2,'   in eleme',
     *       'nt',I6,/)
99989 FORMAT (/' No gradient elements assigned.')
      END
