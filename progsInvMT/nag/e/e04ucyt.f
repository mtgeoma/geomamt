      SUBROUTINE E04UCY(INFORM,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,LDCJ,
     *                  LDCJU,N,NCNLN,CONFUN,OBJFUN,NEEDC,BIGBND,EPSRF,
     *                  CDINT,FDINT,FDCHK,FDNORM,OBJF,XNORM,BL,BU,C,C1,
     *                  CJAC,CJACU,CJDX,DX,GRAD,GRADU,HFORWD,HCNTRL,X,
     *                  WRK1,WRK2,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1092 (JUL 1993).
C     MARK 17 REVISED. IER-1612 (JUN 1995).
C
C     ******************************************************************
C     E04UCY  performs the following...
C     (1)  Computes the objective and constraint values OBJF and C.
C     (2)  Evaluates the user-provided gradients in CJACU and GRADU.
C     (3)  Counts the missing gradients.
C     (4)  Loads the known gradients into GRAD and CJAC.
C     (5)  Checks that the known gradients are programmed correctly.
C     (6)  Computes the missing gradient elements.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written 4-September-1985.
C     This version of E04UCY dated 26-Nov-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CDINT, EPSRF, FDCHK, FDINT, FDNORM,
     *                  OBJF, XNORM
      INTEGER           INFORM, LDCJ, LDCJU, LVLDER, MSGNP, N, NCNLN,
     *                  NFUN, NGRAD, NSTATE
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), CJDX(*), DX(N), GRAD(N),
     *                  GRADU(N), HCNTRL(*), HFORWD(*), USER(*),
     *                  WRK1(N+NCNLN), WRK2(N+NCNLN), X(N)
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
      INTEGER           I, INFOG, INFOJ, J, MODE, NCSET
      LOGICAL           CENTRL, NEEDFD
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, E04UDS, E04XAW, E04XAX, E04XAY, F06DBF,
     *                  F06FBF, F06QFF, F06QHF, X04BAY
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
C     .. Executable Statements ..
C
      INFOG = 0
      INFOJ = 0
      NFDIFF = 0
      NCDIFF = 0
      NCSET = N*NCNLN
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
         IF (MODE.LT.0) GO TO 80
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
      IF (NEEDFD) CALL F06FBF(N,RDUMMY,GRADU,1)
C
      MODE = 2
      CALL OBJFUN(MODE,N,X,OBJF,GRADU,NSTATE,IUSER,USER)
      IF (MODE.LT.0) GO TO 80
C
      CALL DCOPY(N,GRADU,1,GRAD,1)
C
      IF (NEEDFD) THEN
C
C        Count the number of missing gradient elements.
C
         DO 60 J = 1, N
            IF (GRADU(J).EQ.RDUMMY) NFDIFF = NFDIFF + 1
   60    CONTINUE
C
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
                  WRITE (REC,FMT=99996) N - NFDIFF, N, NFDIFF
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
            IF (MODE.LT.0) GO TO 80
            INFOJ = MODE
         END IF
C
         IF (NFDIFF.LT.N) THEN
            CALL E04XAX(MODE,MSGNP,N,BIGBND,EPSRF,EPSPT3,FDCHK,OBJF,
     *                  XNORM,OBJFUN,BL,BU,GRAD,GRADU,DX,X,WRK1,IUSER,
     *                  USER)
            IF (MODE.LT.0) GO TO 80
            INFOG = MODE
         END IF
      END IF
C
      NEEDFD = NCDIFF .GT. 0 .OR. NFDIFF .GT. 0
      IF (NEEDFD) THEN
C        ===============================================================
C        Compute the missing gradient elements.
C        ===============================================================
         CALL E04XAY(MODE,MSGNP,LVLDER,N,NCNLN,LDCJ,LDCJU,BIGBND,EPSRF,
     *               FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,C1,CJDX,
     *               CJAC,CJACU,GRAD,GRADU,HFORWD,HCNTRL,X,DX,IUSER,
     *               USER)
C
         IF (MODE.LT.0) GO TO 80
C
         IF (LFDSET.GT.0) THEN
            CENTRL = LVLDIF .EQ. 2
            CALL E04UDS(CENTRL,MODE,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                  FDINT,FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,
     *                  C1,CJDX,CJAC,CJACU,GRAD,GRADU,HFORWD,HCNTRL,X,
     *                  IUSER,USER)
C
            IF (MODE.LT.0) GO TO 80
         END IF
      END IF
C
      INFORM = INFOJ + INFOG
      RETURN
C
C     The user requested termination.
C
   80 INFORM = MODE
      RETURN
C
C
C     End of  E04UCY. (NPCHKD)
C
99999 FORMAT (/' All Jacobian elements have been set.',/' Derivative l',
     *       'evel increased to ',I4)
99998 FORMAT (/' The user sets ',I6,'   out of',I6,'   Jacobian elemen',
     *       'ts.',/' Each iteration, ',I6,'   Jacobian elements will ',
     *       'be estimated numerically.')
99997 FORMAT (/' All objective gradient elements have been set.',/' De',
     *       'rivative level increased to ',I4)
99996 FORMAT (/' The user sets ',I6,'   out of',I6,'   objective gradi',
     *       'ent elements.',/' Each iteration, ',I6,'   gradient elem',
     *       'ents will be estimated numerically.')
      END
