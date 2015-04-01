      SUBROUTINE E04UPY(INFORM,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,LDCJ,
     *                  LDCJU,LDFJ,LDFJU,M,N,NCNLN,CONFUN,OBJFUN,NEEDC,
     *                  BIGBND,EPSRF,CDINT,FDINT,FDCHK,FDNORM,XNORM,BL,
     *                  BU,C,C1,CJAC,CJACU,CJDX,F,F1,FJAC,FJACU,FJDX,DX,
     *                  HFORWD,HCNTRL,X,WRK1,WRK2,WRK4,IUSER,USER)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1104 (JUL 1993).
C     MARK 17 REVISED. IER-1624 (JUN 1995).
C
C     ******************************************************************
C     E04UPY  performs the following...
C     (1)  Computes the objective and constraint values F and C.
C     (2)  Evaluates the user-provided gradients in CJACU and FJACU.
C     (3)  Counts the missing gradients.
C     (4)  Loads the known gradients into CJAC and FJAC.
C     (5)  Checks that the known gradients are programmed correctly.
C     (6)  Computes the missing gradient elements.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version based on NPCHKD, written 4-September-1985.
C     This version of E04UPY dated 11-May-1988.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CDINT, EPSRF, FDCHK, FDINT, FDNORM,
     *                  XNORM
      INTEGER           INFORM, LDCJ, LDCJU, LDFJ, LDFJU, LVLDER, M,
     *                  MSGNP, N, NCNLN, NFUN, NGRAD, NSTATE
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), CJDX(*), DX(N), F(M), F1(M),
     *                  FJAC(LDFJ,*), FJACU(LDFJU,*), FJDX(M),
     *                  HCNTRL(*), HFORWD(*), USER(*), WRK1(N), WRK2(*),
     *                  WRK4(M), X(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LFDSET, LINES1, LINES2, LVLDIF,
     *                  LVRFYC, NCDIFF, NFDIFF, NOUT
C     .. Arrays in Common ..
      INTEGER           JVERFY(4)
C     .. Local Scalars ..
      INTEGER           I, INFOCJ, INFOFJ, J, MODE, NCSET, NFSET
      LOGICAL           CENTRL, NEEDFD
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          E04UPQ, E04UPR, E04UPW, E04XAW, F06DBF, F06QFF,
     *                  F06QHF, X04BAY
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
C     .. Executable Statements ..
C
      INFOCJ = 0
      INFOFJ = 0
      NFDIFF = 0
      NCDIFF = 0
      NCSET = N*NCNLN
      NFSET = N*M
C
      IF (NCNLN.GT.0) THEN
C        ===============================================================
C        Compute the constraints and Jacobian matrix.
C        ===============================================================
C        If some derivatives are missing, load the Jacobian with dummy
C        values.  Any elements left unaltered after the call to CONFUN
C        must be estimated.  A record of the missing Jacobian elements
C        is stored in  CJACU.
C
         NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 1
C
         IF (NEEDFD) CALL F06QHF('General',NCNLN,N,RDUMMY,RDUMMY,CJACU,
     *                           LDCJU)
C
         CALL F06DBF(NCNLN,(1),NEEDC,1)
C
         MODE = 2
         CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C,CJACU,NSTATE,IUSER,
     *               USER)
         IF (MODE.LT.0) GO TO 100
C
         CALL F06QFF('General',NCNLN,N,CJACU,LDCJU,CJAC,LDCJ)
C
         IF (NEEDFD) THEN
C
C           Count the number of missing Jacobian elements.
C
            DO 40 J = 1, N
               DO 20 I = 1, NCNLN
                  IF (CJACU(I,J).EQ.RDUMMY) NCDIFF = NCDIFF + 1
   20          CONTINUE
   40       CONTINUE
C
            NCSET = NCSET - NCDIFF
            IF (NSTATE.EQ.1) THEN
               IF (NCDIFF.EQ.0) THEN
                  IF (LVLDER.EQ.0) LVLDER = 2
                  IF (LVLDER.EQ.1) LVLDER = 3
                  IF (MSGNP.GT.0) THEN
                     WRITE (REC,FMT=99999) LVLDER
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               ELSE
                  IF (MSGNP.GT.0) THEN
                     WRITE (REC,FMT=99998) NCSET, N*NCNLN, NCDIFF
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     ==================================================================
C     Repeat the procedure above for the objective function.
C     ==================================================================
      NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2
C
      IF (NEEDFD) CALL F06QHF('General',M,N,RDUMMY,RDUMMY,FJACU,LDFJU)
C
      MODE = 2
      CALL OBJFUN(MODE,M,N,LDFJU,X,F,FJACU,NSTATE,IUSER,USER)
      IF (MODE.LT.0) GO TO 100
C
      CALL F06QFF('General',M,N,FJACU,LDFJU,FJAC,LDFJ)
C
      IF (NEEDFD) THEN
         DO 80 J = 1, N
            DO 60 I = 1, M
               IF (FJACU(I,J).EQ.RDUMMY) NFDIFF = NFDIFF + 1
   60       CONTINUE
   80    CONTINUE
C
         NFSET = NFSET - NFDIFF
         IF (NSTATE.EQ.1) THEN
            IF (NFDIFF.EQ.0) THEN
               IF (LVLDER.EQ.0) LVLDER = 1
               IF (LVLDER.EQ.2) LVLDER = 3
               IF (MSGNP.GT.0) THEN
                  WRITE (REC,FMT=99997) LVLDER
                  CALL X04BAY(IPRINT,3,REC)
               END IF
            ELSE
               IF (MSGNP.GT.0) THEN
                  WRITE (REC,FMT=99996) NFSET, N*M, NFDIFF
                  CALL X04BAY(IPRINT,3,REC)
               END IF
            END IF
         END IF
      END IF
C
      NFUN = NFUN + 1
      NGRAD = NGRAD + 1
C
C     ==================================================================
C     Check whatever gradient elements have been provided.
C     ==================================================================
      IF (LVRFYC.GE.0) THEN
         IF (NCSET.GT.0) THEN
            CALL E04XAW(MODE,LVLDER,MSGNP,NCSET,N,NCNLN,LDCJ,LDCJU,
     *                  BIGBND,EPSRF,EPSPT3,FDCHK,XNORM,CONFUN,NEEDC,BL,
     *                  BU,C,C1,CJAC,CJACU,CJDX,DX,WRK2,X,WRK1,IUSER,
     *                  USER)
            IF (MODE.LT.0) GO TO 100
            INFOCJ = MODE
         END IF
C
         IF (NFSET.GT.0) THEN
            CALL E04UPQ(MODE,LVLDER,MSGNP,NFSET,M,N,LDFJ,LDFJU,BIGBND,
     *                  EPSRF,EPSPT3,FDCHK,XNORM,OBJFUN,BL,BU,F,F1,FJAC,
     *                  FJACU,FJDX,DX,WRK4,X,WRK1,IUSER,USER)
            IF (MODE.LT.0) GO TO 100
            INFOFJ = MODE
         END IF
      END IF
C
      NEEDFD = NCDIFF .GT. 0 .OR. NFDIFF .GT. 0
      IF (NEEDFD) THEN
C        ===============================================================
C        Compute the missing gradient elements.
C        ===============================================================
         CALL E04UPR(MODE,MSGNP,LVLDER,M,N,NCNLN,LDCJ,LDCJU,LDFJ,LDFJU,
     *               BIGBND,EPSRF,FDNORM,CONFUN,OBJFUN,NEEDC,BL,BU,C,C1,
     *               CJDX,CJAC,CJACU,F,F1,FJDX,FJAC,FJACU,HFORWD,HCNTRL,
     *               X,DX,IUSER,USER)
C
         IF (MODE.LT.0) GO TO 100
C
         IF (LFDSET.GT.0) THEN
            CENTRL = .FALSE.
            CALL E04UPW(CENTRL,MODE,LDCJ,LDCJU,LDFJ,LDFJU,M,N,NCNLN,
     *                  BIGBND,CDINT,FDINT,FDNORM,CONFUN,OBJFUN,NEEDC,
     *                  BL,BU,C,C1,CJDX,CJAC,CJACU,F,F1,FJDX,FJAC,FJACU,
     *                  HFORWD,HCNTRL,X,IUSER,USER)
C
            IF (MODE.LT.0) GO TO 100
         END IF
      END IF
C
      INFORM = INFOCJ + INFOFJ
C
      RETURN
C
C     The user requested termination.
C
  100 INFORM = MODE
      RETURN
C
C
C     End of  E04UPY.  (NLCHKD)
C
99999 FORMAT (/' All constraint Jacobian elements have been set.',/' D',
     *       'erivative level increased to ',I4)
99998 FORMAT (/' The user sets ',I6,'   out of',I6,'   constraint Jaco',
     *       'bian elements.',/' Each iteration, ',I6,'   constraint J',
     *       'acobian elements will be estimated.')
99997 FORMAT (/' All objective Jacobian elements have been set.',/' De',
     *       'rivative level increased to ',I4)
99996 FORMAT (/' The user sets ',I6,'   out of',I6,'   objective Jacob',
     *       'ian elements.',/' Each iteration, ',I6,'   objective Jac',
     *       'obian elements will be estimated.')
      END
