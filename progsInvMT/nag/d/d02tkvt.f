      SUBROUTINE D02TKV(MSING,XI,XIOLD,Z,DMZ,DELZ,DELDMZ,G,W,V,RHS,DMZO,
     *                  INTEGS,IPVTG,IPVTW,RNORM,MODE,FFUN,FJAC,GAFUN,
     *                  GAJAC,GBFUN,GBJAC,GUESS,IGUESS,M,MAXORD,NEQ,DGR,
     *                  LDGR,DGZ,DF,F,KD,MSTAR,NLBC,SETRES,KCOL,RHO,
     *                  COEF,B,ACOL,N,NOLD,NMAX,ZVAL,ZVAL1)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C*********************************************************************
C
C   Purpose
C      This routine controls the set up and solution of a linear
C      system of collocation equations.
C      The matrix  G  is cast into an almost block diagonal
C      form by an appropriate ordering of the columns and solved
C      using NAG Fortran Library routines F01LHF and F04LHF.
C      The matrix is composed of N+2 blocks. The i-th block has size
C         INTEGS(1,i) rows   by   INTEGS(2,i) columns.
C      For i=2,3,...,N+1, the i-th block contains the rows the
C      linearized collocation equations for the (i-1)th interval,
C      condensed as described in [2]. Hence its size is
C      (MSTAR x 2*MSTAR). The first block contains the linearized
C      left hand boundary conditions (hence of size NLBC x MSTAR)
C      and the last (N+2)-th block contains the linearized
C      right hand boundary conditions (hence of size
C      (MSTAR-NLBC) x MSTAR). INTEGS(3,i) contains the number of
C      columns of overlap between the i-th and (i+1)-th blocks,
C      which is always MSTAR in this application.
C      The right hand side vector is put into RHS and the solution
C      is returned in DELZ and DELDMZ.
C
C      It operates operates according to one of 5 modes:
C      MODE = 0 - set up the collocation matrices V , W , G
C                 and the right hand side RHS, and solve.
C                 (For linear problems only.)
C      MODE = 1 - set up the collocation matrices V , W , G
C                 and the right hand sides RHS and DMZO,
C                 and solve. Also set up INTEGS.
C                 (First iteration of nonlinear problems only).
C      MODE = 2 - set up RHS only and compute its norm.
C      MODE = 3 - set up V, W, G only and solve system.
C      MODE = 4 - perform forward and backward substitution only
C                 (do not set up the matrices nor form the RHS).
C
C   Arguments
C      MSING  - indicates if linear system was singular
C      XI     - current mesh
C      XIOLD  - previous mesh
C      Z      - solution values on each subinterval and boundary
C               conditions
C      DMZ    - high order solution derivatives at each collocation
C               point of each interval for each equation
C      DELZ   - solution of linear system (hence, correction to Z)
C      DELDMZ - uncondensed solution (hence, correction ot DMZ)
C      G      - condensed matrix and its decomposition
C      W      - submatrices used in condensation for each interval
C      V      - submatrices used in condensation for each interval
C      RHS    - actual right hand side
C      DMZO   -
C      INTEGS - array to describe structure of linear system in G
C      IPVTG  - pivot index for decompostion of G
C      IPVTW  - pivot index for each submatrix W
C      RNORM  - norm of RHS
C      MODE   - mode of operation for LSYSLV
C      FFUN   - procedure to evaluate ODEs
C      FJAC  - procedure to evaluate Jacobian of ODEs
C      GSUB   - procedure to evaluate boundary conditions
C      DGSUB  - procedure to evaluate Jacobian of boundary conditions
C      GUESS  - procedure to evaluate initial solution (if supplied)
C      IGUESS - mode of provision of initial solution (1 via GUESS, or 0
C      M      - orders of ODEs
C      MAXORD - maximal order of ODEs
C      NEQ    - number of ODEs
C      DGR    - array for Jacobian of boundary conditions
C      LDGR   - leading dim of DGR
C      DGZ    - inner product of Jacobian and solution values
C      DF     - array for Jacobian of ODEs
C      F      - array for ODEs
C      KD     - dimension of workspaces (NEQ*KCOL)
C      MSTAR  - number of solution components (sum(M(i),i=1,NEQ))
C      NLBC   - number of boundary conditions at left hand side
C      SETRES - indicates if DMZO (initial residual) to be evalauted
C      KCOL   - number of collocation points
C      RHO    - collocation points
C      COEF   - mesh dependent monomial coefficients
C      B      - rk-basis coefficients for each interval
C      ACOL   - rk-basis coefficients for collocation points
C      N      - number of points in current mesh
C      NOLD   - number of points in previous mesh
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine LSYSLV)
C
C*********************************************************************
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RNORM
      INTEGER           IGUESS, KCOL, KD, LDGR, MAXORD, MODE, MSING,
     *                  MSTAR, N, NEQ, NLBC, NMAX, NOLD
      LOGICAL           SETRES
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOL(28,7), B(28), COEF(KCOL,KCOL),
     *                  DELDMZ(N*KCOL*NEQ), DELZ(MSTAR*(N+1)),
     *                  DF(NEQ,NEQ,MAXORD), DGR(LDGR,NEQ,MAXORD),
     *                  DGZ(MSTAR), DMZ(NMAX*KCOL*NEQ),
     *                  DMZO(N*KCOL*NEQ), F(NEQ),
     *                  G((N*2+1)*MSTAR*MSTAR), RHO(KCOL),
     *                  RHS(N*KCOL*NEQ+MSTAR), V(KD,MSTAR,N),
     *                  W(KD,KD,N), XI(N+1), XIOLD(NOLD+1),
     *                  Z(MSTAR*(NMAX+1)), ZVAL(MSTAR),
     *                  ZVAL1(NEQ,MAXORD)
      INTEGER           INTEGS(3,N+2), IPVTG(MSTAR*(N+1)), IPVTW(KD,N),
     *                  M(NEQ)
C     .. Subroutine Arguments ..
      EXTERNAL          FFUN, FJAC, GAFUN, GAJAC, GBFUN, GBJAC, GUESS
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HRHO, VALUE, XCOL, XII
      INTEGER           I, I1, IDMZ, IDMZO, IG, INDEX, IOLD, IRHS, IZ,
     *                  J, J1, JCOL, JJ, L, LENG, LGIBLK, NCOL, NDMZ,
     *                  NRBC, NZ
      LOGICAL           GENRHS, GENROW, GENZVL
C     .. Local Arrays ..
      DOUBLE PRECISION  AT(28), DMVAL(20), DUMMY(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           D02TKS
      EXTERNAL          X02AJF, D02TKS
C     .. External Subroutines ..
      EXTERNAL          D02TKM, D02TKN, D02TKP, D02TKQ, D02TKR, D02TKT,
     *                  F01LHF, F04LHF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
      NZ = MSTAR*(N+1)
      NDMZ = KCOL*NEQ*N
      LENG = (N*2+1)*MSTAR*MSTAR
      IF (MODE.EQ.4) GO TO 340
C
C Linear problem initialization
C
      IF (MODE.EQ.0) THEN
         DO 20 I = 1, MSTAR
            ZVAL(I) = 0.D0
   20    CONTINUE
      END IF
C
C Initialization
C
      IDMZ = 1
      IDMZO = 1
      IRHS = 1
      IG = 1
      IOLD = 1
      NCOL = 2*MSTAR
      RNORM = 0.D0
C
C Build integs (describing block structure of matrix)
C
      NRBC = MSTAR - NLBC
      INTEGS(1,1) = NLBC
      INTEGS(2,1) = MSTAR
      INTEGS(3,1) = MSTAR
      INTEGS(1,N+2) = NRBC
      INTEGS(2,N+2) = MSTAR
      DO 40 I = 2, N + 1
         INTEGS(1,I) = MSTAR
         INTEGS(2,I) = NCOL
         INTEGS(3,I) = MSTAR
   40 CONTINUE
      LGIBLK = NCOL*MSTAR
C
C Zero the matrices to be computed if necessary
C
      IF (MODE.NE.2) THEN
         DO 100 L = 1, N
            DO 80 J = 1, KD
               DO 60 I = 1, KD
                  W(I,J,L) = 0.D0
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
      END IF
C
C Set up the linear equations block by block
C
      GENZVL = MODE .EQ. 1 .OR. MODE .EQ. 2 .OR. MODE .EQ. 3
      GENRHS = MODE .EQ. 0 .OR. MODE .EQ. 1 .OR. MODE .EQ. 2
      GENROW = MODE .EQ. 0 .OR. MODE .EQ. 1 .OR. MODE .EQ. 3
C
C Construct first block and/or corresponding piece of RHS
C    due to left hand boundary conditions
C
C      SETRES = NONLIN .NE. 0 .AND. ITER .EQ. 0
      IF (GENZVL) THEN
         IF (IGUESS.EQ.1) THEN
            CALL GUESS(XI(1),NEQ,M,ZVAL1,DMVAL)
C            JCOL = 0
C            DO 160 I = 1, NEQ
C               DO 140 J = 1, M(I)
C                  JCOL = JCOL + 1
C                  ZVAL1(I,J) = ZVAL(JCOL)
C  140          CONTINUE
C  160       CONTINUE
         ELSE
C
C Note, regardless of mode (1, >1) left hand end is always positioned
C at the beginning of Z
C
            JCOL = 0
            DO 140 I = 1, NEQ
               DO 120 J = 1, M(I)
                  JCOL = JCOL + 1
                  ZVAL1(I,J) = Z(JCOL)
  120          CONTINUE
  140       CONTINUE
         END IF
      END IF
C
      IF (GENRHS) THEN
         CALL GAFUN(ZVAL1,NEQ,M,NLBC,RHS(NDMZ+1))
         DO 160 I = 1, NLBC
            RNORM = RNORM + RHS(NDMZ+I)**2
            RHS(NDMZ+I) = -RHS(NDMZ+I)
  160    CONTINUE
      END IF
      IF (GENROW) THEN
         CALL D02TKN(G(1),NLBC,ZVAL1,DGZ,GAJAC,M,MAXORD,NEQ,MSTAR,DGR,
     *               SETRES)
      END IF
C
C Construct last block and/or corresponding piece of RHS
C    due to right hand boundary conditions
C
      IF (GENZVL) THEN
         IF (IGUESS.EQ.1) THEN
            CALL GUESS(XI(N+1),NEQ,M,ZVAL1,DMVAL)
C            JCOL = 0
C            DO 260 I = 1, NEQ
C               DO 240 J = 1, M(I)
C                  JCOL = JCOL + 1
C                  ZVAL1(I,J) = ZVAL(JCOL)
C  240          CONTINUE
C  260       CONTINUE
         ELSE
            IF (MODE.EQ.1) THEN
               JCOL = NOLD*MSTAR
            ELSE
               JCOL = N*MSTAR
            END IF
            DO 200 I = 1, NEQ
               DO 180 J = 1, M(I)
                  JCOL = JCOL + 1
                  ZVAL1(I,J) = Z(JCOL)
  180          CONTINUE
  200       CONTINUE
         END IF
      END IF
C
      IF (GENRHS) THEN
         CALL GBFUN(ZVAL1,NEQ,M,NRBC,RHS(NDMZ+NLBC+1))
         DO 220 I = NLBC + 1, NLBC + NRBC
            RNORM = RNORM + RHS(NDMZ+I)**2
            RHS(NDMZ+I) = -RHS(NDMZ+I)
  220    CONTINUE
      END IF
      IF (GENROW) THEN
         IG = 1 + NLBC*MSTAR + N*LGIBLK
         CALL D02TKN(G(IG),NRBC,ZVAL1,DGZ(NLBC+1),GBJAC,M,MAXORD,NEQ,
     *               MSTAR,DGR,SETRES)
      END IF
C
C Construct intermediate blocks and/or corresponding pieces of RHS
C    due to collocation subintervals
C
      IG = 1 + NLBC*MSTAR
      DO 320 I = 1, N
         XII = XI(I)
         H = XI(I+1) - XI(I)
C
C Assemble collocation equations
C
         DO 300 J = 1, KCOL
            HRHO = H*RHO(J)
            XCOL = XII + HRHO
C
C XCOL corresponds to a collocation (interior) point.
C    Build the corresponding  neq  equations.
C
            IF (GENZVL) THEN
               IF (IGUESS.EQ.1) THEN
                  CALL GUESS(XCOL,NEQ,M,ZVAL1,DMZO(IRHS))
               ELSE
                  IF (MODE.EQ.1) THEN
                     IOLD = D02TKS(XCOL,XIOLD,NOLD)
                     IZ = (IOLD-1)*MSTAR + 1
                     IDMZ = (IOLD-1)*KD + 1
                     CALL D02TKT(XCOL,ZVAL,AT,COEF,XIOLD(IOLD),
     *                           XIOLD(IOLD+1),Z(IZ),DMZ(IDMZ),KCOL,NEQ,
     *                           MAXORD,M,MSTAR,.TRUE.,.TRUE.,DMZO(IRHS)
     *                           )
                  ELSE
                     IZ = (I-1)*MSTAR + 1
                     IDMZ = (I-1)*KD + 1
                     CALL D02TKT(XCOL,ZVAL,ACOL(1,J),COEF,XI(I),XI(I+1),
     *                           Z(IZ),DMZ(IDMZ),KCOL,NEQ,MAXORD,M,
     *                           MSTAR,.FALSE.,.FALSE.,DUMMY)
                  END IF
C
                  JCOL = 0
                  DO 260 I1 = 1, NEQ
                     DO 240 J1 = 1, M(I1)
                        JCOL = JCOL + 1
                        ZVAL1(I1,J1) = ZVAL(JCOL)
  240                CONTINUE
  260             CONTINUE
               END IF
            END IF
C
            IF (GENRHS) THEN
               IF (MODE.EQ.0) THEN
                  CALL FFUN(XCOL,ZVAL1,NEQ,M,RHS(IRHS))
                  IRHS = IRHS + NEQ
               ELSE
                  CALL FFUN(XCOL,ZVAL1,NEQ,M,F)
                  DO 280 JJ = 1, NEQ
                     IF (MODE.EQ.1) THEN
                        VALUE = DMZO(IRHS) - F(JJ)
                     ELSE
                        VALUE = DMZ(IRHS) - F(JJ)
                     END IF
                     RHS(IRHS) = -VALUE
                     RNORM = RNORM + VALUE**2
                     IRHS = IRHS + 1
  280             CONTINUE
               END IF
            END IF
            IF (GENROW) THEN
C
C Do condensation
C
               CALL D02TKM(XCOL,HRHO,J,KCOL,W(1,1,I),V(1,1,I),IPVTW(1,I)
     *                     ,KD,ZVAL1,DF,ACOL(1,J),DMZO(IDMZO),FJAC,
     *                     MSING,M,NEQ,MAXORD,MSTAR,SETRES)
               IF (MSING.NE.0) RETURN
            END IF
  300    CONTINUE
C
C Build global bvp matrix G
C
         IF (GENROW) CALL D02TKP(H,B,G(IG),MSTAR,V(1,1,I),M,MAXORD,NEQ,
     *                           KCOL)
         IG = IG + LGIBLK
         IDMZ = IDMZ + KD
         IF (MODE.EQ.1) IDMZO = IDMZO + KD
  320 CONTINUE
C
C Assembly process completed
C
      IF (MODE.EQ.1 .OR. MODE.EQ.2) THEN
         RNORM = SQRT(RNORM/DBLE(NZ+NDMZ))
         IF (MODE.EQ.2) RETURN
      END IF
C
C Solve the linear system.
C
C Matrix decomposition
C
      MSING = 1
      CALL F01LHF(NZ,N+2,INTEGS,G,LENG,IPVTG,1.0D-3*X02AJF(),INDEX,
     *            MSING)
C
C Check for singular matrix
C
      IF (MSING.NE.0) THEN
         MSING = -1
         RETURN
      END IF
C
C Perform forward and backward substitution for mode=4 only.
C
  340 CONTINUE
      DO 360 L = 1, NDMZ
         DELDMZ(L) = RHS(L)
  360 CONTINUE
C
      NRBC = MSTAR - NLBC
C
C Account for condensation
C
      DO 380 I = 1, NLBC
         DELZ(I) = RHS(NDMZ+I)
  380 CONTINUE
      IZ = NLBC + 1
      IDMZ = 1
      DO 400 I = 1, N
         H = XI(I+1) - XI(I)
         CALL D02TKQ(H,B,W(1,1,I),IPVTW(1,I),KD,DELZ(IZ),DELDMZ(IDMZ),M,
     *               MAXORD,NEQ,KCOL,MSTAR)
         IZ = IZ + MSTAR
         IDMZ = IDMZ + KD
  400 CONTINUE
      DO 420 I = NLBC + 1, NLBC + NRBC
         DELZ(N*MSTAR+I) = RHS(NDMZ+I)
  420 CONTINUE
C
C Perform forward and backward substitution for mode=0,2, or 3.
C
      CALL F04LHF('n',NZ,N+2,INTEGS,G,LENG,IPVTG,DELZ,NZ,1,MSING)
C
C Finaly find deldmz
C
      CALL D02TKR(KD,MSTAR,N,V,DELZ,DELDMZ)
C
      IF (MODE.NE.1) RETURN
C
      DO 440 L = 1, NDMZ
         DMZ(L) = DMZO(L)
  440 CONTINUE
C
C Account for condensation
C
      DO 460 I = 1, NLBC
         Z(I) = DGZ(I)
  460 CONTINUE
      IZ = NLBC + 1
      IDMZ = 1
      DO 480 I = 1, N
         H = XI(I+1) - XI(I)
         CALL D02TKQ(H,B,W(1,1,I),IPVTW(1,I),KD,Z(IZ),DMZ(IDMZ),M,
     *               MAXORD,NEQ,KCOL,MSTAR)
         IZ = IZ + MSTAR
         IDMZ = IDMZ + KD
  480 CONTINUE
      DO 500 I = NLBC + 1, NLBC + NRBC
         Z(N*MSTAR+I) = DGZ(I)
  500 CONTINUE
C
      CALL F04LHF('n',NZ,N+2,INTEGS,G,LENG,IPVTG,Z,NZ,1,MSING)
C
C Finally find dmz
C
      CALL D02TKR(KD,MSTAR,N,V,Z,DMZ)
C
      RETURN
      END
