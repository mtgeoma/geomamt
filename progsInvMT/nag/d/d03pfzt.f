      SUBROUTINE D03PFZ(PDEFCN,BNDARY,DUMFCN,DUMBND,NEQ,NEQMAX,T,TOUT,Y,
     *                  RTOL,ATOL,ITOL,ITRACE,SNORM,MATZ,RESID,MONITR,
     *                  WKMON,NWKMON,WD,NWD,WKJAC,LENWJ,IW,NIW,OPT,
     *                  LDERIV,ITASK,WKRES,NWKRES,PDEFN,FLXPFF,FLXPLF,
     *                  BNDY,ODEFN,MONITF,IND,IFAIL1,IA,NIA,JA,NJA,
     *                  JCEVL,JACFUL,JACBND,JACSPS,IRES,IXFIX,NIXFIX,NW)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C----------------------------------------------------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C          PARAMETER LIST
C          --------------
C  The following variables are used to communicate with D03PFZ.
C  These variables should be declared in the calling program as follows
C
C      INTEGER  ITRACE, NWD, NIW, IW(NIW), ITOL, ITASK, M, IA, JA,
C    1          IFAIL1, NEQ, NEQMAX, NWKRES, IND, IA(NIA), JA(NJA)
C      DOUBLE PRECISION T, TOUT, Y(NEQMAX), WD(NWD), RTOL(*),
C    1          ATOL(*), OPT(30)
C      CHARACTER*1 SNORM, MATZ, METOD1, METOD2, JCEVAL
C      LOGICAL LBLOCK, LDERIV(2)
C      EXTERNAL RESID, MONITR, JAC, PDEFN, FLXPFF,FLXPLF, BNDY, ODEFN,
C               MONITF
C
C  The values that should be assigned to these variables are
C
C  NEQ - INTEGER
C    Number of equations
C
C  NEQMAX - INTEGER
C    Upper bound on number of number of equations
C
C  T - DOUBLE PRECISION
C    Initial value of time (the independent variable) (input). On output
C    it contains the time at which a computed solution is returned.
C
C  TOUT - DOUBLE PRECISION
C    Output time
C
C  Y - DOUBLE PRECISION
C    Array of length (NEQMAX)
C    On input it must contain the initial values of the neqmax variables
C    On output it contains the computed solution at the final output
C    time TOUT. If the solution fails before TOUT, or is
C    stopped by the user in subroutine RESID or MONITR (see below),
C    Y contains the last-computed solution, and T contains the time that
C    was reached .
C
C  RTOL - DOUBLE PRECISION
C   Array of length (*) containing the relative error tolerances for
C   each variable (input). The error bounds used in controlling the
C   error are calculated as
C
C              EWT(I) = RTOL(I) * ABS(Y(I)) + ATOL(I)
C
C   Frequently the same value of rtol can be used for all variables a
C   a value of 1d-4 should give moderate accuracy.
C   To check validity of results do a repeat run with different RTOL
C   values, and see how much the solution changes.
C
C  ATOL - DOUBLE PRECISION
C   Array of length (*) containing the absolute error tolerances for
C   variable (input) (see RTOL).
C   Different values are required for different variables
C   if there are large differences in scaling. For each variable,
C   a value of atol should be chosen that is a few orders
C   of magnitude smaller than a typical value for the variable.
C
C  ITOL - INTEGER
C   Controls form of local error test.
C       = 1 - RTOL, ATOL both scalar
C       = 2 - RTOL: VECTOR, ATOL: scalar
C       = 3 - RTOL: SCALAR, ATOL: vector
C       = 4 - RTOL, ATOL both vector
C
C  ITRACE - INTEGER
C   Specifies level of trace information required from D02NAG, permitted
C   values  0,1,2 or 3 (input). If problems occur
C   turn on the trace to obtain diagnostic information.
C   the output is automatically written to the nag trace channel which
C   is assigned to channel whose number is contained in idev by the call
C       CALL X04ABF( 1, IDEV)
C
C  SNORM
C   Character variable specifying type of norm to be used by d03nag
C   in controlling the local error (input) .
C   possible values are  '1' and '2' (L1 and L2 norm respectively).
C
C  MATZ
C   Character variable specifying type of matrix algebra
C   required (input). The possible choices are
C      'F', 'B', 'S'.
C   Use 'F' for small problems, and unless
C   the system is known to have a banded or sparse coupling pattern.
C
C
C   RESID
C   Name of a subroutine supplied by the user to define the equations.
C   it must be declared as external in the user's calling program.
C   the user can stop the integration from resid by setting its
C   parameter ires to a particular value.
C
C   MONITR
C   Name of a subroutine needed if the user wishes to provide output
C   at each time step, to perform intermediate calculation or to control
C   (stop) the integration.
C   It must be declared as external in the user's calling program.
C   If the user has no requirements for output/action
C   At each time step s/he may select a 'dummy' routine.
C
C  WKMON
C   DOUBLE PRECISION array used as workspace for the monitr routine.
C   Assumed to be dimensioned as WKMON(NWKMON)
C
C  NWKMON
C   Size of the workspace required by the monitr routine must be an
C   integer greater than or equal to zero.
C
C  WD
C   double precision workspace of length NWD.
C
C  NWD
C   (note use of blend is not allowed by this code)
C   Size of W (INPUT). The required size of workspace is
C                      NWD  >=.  (11*NEQMAX+50)
C
C    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    * Note in the case when the NDASSL integrator is used the size *
C    *  of the workspace must be increased by 2*NEQMAX              *
C    *  the size of the integrator workspace can be reduced if the  *
C    *  THETA method is selected (OPT(1) = 2) by   2*NEQMAX   or if *
C    *  the BDF integrator is used with a maximum order of less     *
C    *  than five by (5-OPT(2))*NEQMAX                              *
C    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  WKJAC(LENWJ)
C   DOUBLE PRECISION workspace used by the linear algebra routines
C
C  LENWJ
C   The size of WKJAC depends on the type of matrix algebra selected :
C
C           FULL   : LENWJ >=.  NEQMAX*NEQMAX + 2
C
C           BANDED : LENWJ >=.  (2*ML+MU+1)*NEQMAX + 2
C
C                      ML and MU are lower and upper bandwidths, set in
C                      the setup routine for the space discretisation.
C
C           SPARSE : LENWJ >=.  4*NEQMAX + 11*NEQMAX/2
C                   This is an initial estimate - the user may need
C                   to provide more storage in which case a message is
C                   written to the terminal.
C
C  IW
C   INTEGER workspace of length NIW.
C
C  NIW
C   Size of IW (INPUT). The required size for matz = 'F'
C                        is          NIW =             24
C                       or in the 'B' case:
C                                    NIW =  2*NEQMAX + 23
C                       or in the 'S' case:
C                                    NIW =  25*NEQMAX + 23
C
C  OPT - DOUBLE PRECISION array of size 30
C   Contains various optional and other inputs:
C
C      OPT           *     BDF         *       THETB2     *  NDASSL
C  (array position)  *  (integrator)   *    (integrator)  * (integrator)
C  --------------------------------------------------------------------
C      (1)           *     1.D0        *       2.D0       *   3.D0
C      (2)           *     5.D0        *  0.51D0 - 0.99D0 *   5.D0
C      (3)           *   1.D0(NEWTON)  *    1.D0(NEWTON)  *
C                    *   2.D0(F/ITER)  *    2.D0(F/ITER)  *
C      (4)           *     1.D0        *    3.D0(SWITCH)  *
C                    *     2.D0        *    4.D0(NOSWCH)  *
C      (5)=CONST(1)  *                 *                  *
C       !    !       *      See integrator documentation  *
C       !    !       *                 "                  *
C      (10)=CONST(6) *                 "                  *
C      (11)=TCRIT    *      See setup documentation       *
C      (12)=HMIN     *                 "                  *
C      (13)=HMAX     *                 "                  *
C      (14)=H0       *                 "                  *
C      (15)=MAXSTP   *                 "                  *
C      (16)=MXHNIL   *                 "                  *
C      (17)=SENS     *                 "                  *
C      (18)=U        *                 "                  *
C      (19)=ETA      *                 "                  *
C      (20)=ISPLIT   *                 "                  *
C      (21)=LBLOCK   *                 "                  *
C      (22)=JCEVAL   *                 "                  *
C   ------------------------------------------------------------------
C
C  LDERIV(2) - LOGICAL
C   LDERIV(1)=.TRUE. - user has supplied both an initial Y and
C   YDOT, otherwise only initial Y is assumed to have been supplied.
C   LDERIV(2)=.TRUE. - means a modified newton method should be used to
C   evaluate iniial Y and YDOT.
C   LDERIV(2)=.FALSE. - A functional iteration is used.
C
C  ITASK - INTEGER
C   Specifies the task to be performed.
C        = 1 Normal computation of output values Y(T) at T=TOUT.
C        = 2 One step and return
C        = 3 Stop at first internal integration point at or beyond
C            T=TOUT and return
C        = 4 Normal computation of output values Y(T) at T=TOUT but
C            without overshooting T=TCRIT
C        = 5 One step and return, without passing TCRIT
C
C
C  WKRES(NWKRES) -DOUBLE PRECISION array
C   Used as the workspace by the routine that defines the d.a.e. system
C
C  NWKRES - INTEGER
C   Specifies the dimension of array WKRES.
C
C     where in the case of the D03PLF routine the discretisation
c     workspace size is given by:
C
C      NWKRES  = (17 + 2*NPDE + 6*NXI + MAXNPT) * NPDE +
C                 7 * MAXNPT + NXI + NV + 1 + NXFIX
C
C   Set NWKRES = 1, if this option not required.
C
C  PDEFN
C   Name of subroutine used to define pde being solved . This name
C   is only passed into the RESID routine. PDEFN must be declared
C   as EXTERNAL in the calling program.
C
C  FLXPFF,FLXPLF
C   Names of subroutine used to supply the convective fluxes. (one
C   actual and one dummy). These names are passed only to the
C   RESID routine.
C
C  BNDY
C   Name of subroutine used to define b.c. of pde being solved. This
C   name is only passed into the RESID routine. BNDY must be declared
C   as EXTERNAL in the calling program.
C
C  ODEFN
C   Name of subroutine used to define ode of ode/pde being solved. This
C   name is only passed into the RESID routine. BNDY must be declared
C   as EXTERNAL in the calling program.
C
C   MONITF
C   Name of subroutine used to define the monitor function used in
C   remeshing. This name is only passed into the monitor routine and
C   is not used elsewhere. MONITF must be declared as EXTERNAL in the
C   calling program.
C
C  IND
C   Error/status indicator in the call to D03PFZ.
C
C  IFAIL1 - INTEGER
C   On entry
C           = 0
C       or  = 1
C
C  NIA
C   Size of array IA
C
C  IA(NIA)
C   Array used to pass sparsity pattern information.
C
C  NJA
C   Size of array JA
C
C  JA(NJA)
C   Array used to pass sparsity pattern information.
C
C  JCEVAL - CHARACTER*6
C   Indicates the technique to be used to compute the
C   jacobian.
C               JCEVAL = 'N' - Jacobian evaluated numerically.
C               JCEVAL = 'A' - Jacobian evaluated analitically.
C  Sparse case also    = 'S' - Structural.
C
C  JACFUL
C   The name of the analytic jacobian routine when MATZ = 'F'
C   must be declared as external in the calling program.
C   Use the nag dummy routine D02NGZ otherwise.
C
C  JACBND
C   The name of the analytic jacobian routine when MATZ = 'B'
C   must be declared as external in the calling program.
C   Use the nag dummy routine D02NHZ otherwise .
C
C  JACSPS
C   The name of the analytic jacobian routine when MATZ = 'S'
C   must be declared as external in the calling program.
C   Use the nag dummy routine D02NJZ otherwise .
C
C----------------------------------------------------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     .. Parameters ..
      INTEGER           HU, H, HMIN1, HMAX1, TN, EL0
      PARAMETER         (HU=15,H=16,HMIN1=17,HMAX1=18,TN=19,EL0=20)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TOUT
      INTEGER           IFAIL1, IND, IRES, ITASK, ITOL, ITRACE, LENWJ,
     *                  NEQ, NEQMAX, NIA, NIW, NIXFIX, NJA, NW, NWD,
     *                  NWKMON, NWKRES
      CHARACTER         JCEVL, MATZ, SNORM
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), OPT(30), RTOL(*), WD(NWD),
     *                  WKJAC(LENWJ), WKMON(NWKMON), WKRES(NWKRES),
     *                  Y(NEQMAX)
      INTEGER           IA(NIA), IW(NIW), IXFIX(*), JA(NJA)
      LOGICAL           LDERIV(2)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, BNDY, DUMBND, DUMFCN, FLXPFF, FLXPLF,
     *                  JACBND, JACFUL, JACSPS, MONITF, MONITR, ODEFN,
     *                  PDEFCN, PDEFN, RESID
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           IIFLAG, IOVFLO, LDUM, LENINF, LENJP, LENRW,
     *                  LTRACE, ML, MU, NY2DIM
      LOGICAL           REMESH
C     .. Local Scalars ..
      DOUBLE PRECISION  ETA, H0, HMAX, HMIN, SENS, TCRIT, THETA, U
      INTEGER           I, ICALL, IDEV, IGROW, IIFAIL, IMON, INFPOS,
     *                  INLN, IPLACE, IPOSRW, IPOSYD, IPOSYS, IREVCM,
     *                  ISPLIT, ITOTW, J, JACPOS, JTRACE, LACOR, LENYD,
     *                  LENYS, LEWT, LIWREQ, LIWUSD, LRWREQ, LRWUSD,
     *                  LSAVR, MAXORD, MAXSTP, MXHNIL, NBLOCK, NGP, NLU,
     *                  NNZ
      LOGICAL           LBLOCK, PETZLD
      CHARACTER         JCEVAL, METOD1, METOD2, TSNORM
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  CONST(6)
C     .. External Subroutines ..
      EXTERNAL          D02MWY, D02NNF, D02NNQ, D02NRF, D02NSF, D02NTF,
     *                  D02NUF, D02NVF, D02NXF, E04UDU, X04ABF
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Common blocks ..
      COMMON            /AD02NM/LTRACE, LDUM
      COMMON            /AD03PC/ML, MU
      COMMON            /FD02NM/DUNFLO, UROUND, IOVFLO
      COMMON            /JD03PC/REMESH
      COMMON            /WD03PC/LENRW, NY2DIM, LENINF, LENJP
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /AD03PC/, /FD02NM/, /JD03PC/, /WD03PC/,
     *                  /AD02NM/, /XD03PC/
C     .. Executable Statements ..
C
      CALL X04ABF(0,IDEV)
      IF (IND.EQ.0) THEN
         JCEVAL = JCEVL
         IND = 1
         IF (OPT(1).EQ.1.D0) THEN
            MAXORD = INT(OPT(2))
            NY2DIM = MAXORD + 1
            IF (OPT(3).EQ.1.D0) THEN
               METOD1 = 'N'
            ELSE IF (OPT(3).EQ.2.D0) THEN
               METOD1 = 'F'
            ELSE
               GO TO 220
            END IF
C
            IF (OPT(4).EQ.1.D0) THEN
               PETZLD = .TRUE.
            ELSE IF (OPT(4).EQ.2.D0) THEN
               PETZLD = .FALSE.
            ELSE
               GO TO 240
            END IF
         ELSE IF (OPT(1).EQ.2.D0) THEN
            NY2DIM = 4
            THETA = OPT(5)
            IF (OPT(6).EQ.1.D0) THEN
               METOD1 = 'N'
            ELSE IF (OPT(6).EQ.2.D0) THEN
               METOD1 = 'F'
            ELSE
               GO TO 220
            END IF
C
            IF (OPT(7).EQ.1.D0) THEN
               METOD2 = 'S'
            ELSE IF (OPT(7).EQ.2.D0) THEN
               METOD2 = 'N'
            ELSE
               GO TO 260
            END IF
         END IF
C
         DO 20 I = 1, 6
            CONST(I) = 0.0D0
   20    CONTINUE
C
         TCRIT = OPT(11)
         HMIN = OPT(12)
         HMAX = OPT(13)
         H0 = OPT(14)
C
         MAXSTP = INT(OPT(15))
         MXHNIL = 0
         SENS = 0.D0
         U = OPT(29)
         ETA = OPT(30)
         ISPLIT = 0
C
         LBLOCK = .FALSE.
C
C
C ... Work out the positions of the  arrays in W and IW ...
C
         INFPOS = 1
         LENINF = 23
C
C     ... The JACPVT array ...
C
         JACPOS = 24
         IF (MATZ.EQ.'F') THEN
            LENJP = 1
         ELSE IF (MATZ.EQ.'B') THEN
            LENJP = NEQMAX
         ELSE IF (MATZ.EQ.'S') THEN
            LENJP = NIW - 23
         END IF
C
C    ... The RWORK array:
C
         IPOSRW = 1
         LENRW = 50 + 4*NEQMAX
C
C    ... The YDOTI array:
C
         IPOSYD = IPOSRW + LENRW
         LENYD = NEQMAX
C
C    ... The YSAVE array:
C
         IPOSYS = IPOSYD + LENYD
         LENYS = NY2DIM*NEQMAX
C
C    ... Check the WKJAC array dimension
C
         IF (MATZ.EQ.'F') THEN
            I = NEQMAX*NEQMAX + NEQMAX
         ELSE IF (MATZ.EQ.'B') THEN
            I = (2*ML+MU+1)*NEQMAX
         ELSE IF (MATZ.EQ.'S') THEN
            I = 4*NEQMAX + 11*NEQMAX/2
         END IF
C
         IF (I.GT.LENWJ) THEN
            GO TO 320
         END IF
C
C ... Check that the workspaces are big enough ...
C
         ITOTW = LENYD + LENYS + LENRW
C
C ... Assign spare storage to Jacobian formation ...
C
         IF (ITOTW.GT.NWD) THEN
            GO TO 300
         END IF
C
C ... Set up the linear algebra routines ...
C
         IF (MATZ.EQ.'F') THEN
            CALL D02NSF(NEQ,NEQMAX,JCEVAL,LENWJ,WD(IPOSRW),IFAIL1)
         ELSE IF (MATZ.EQ.'B') THEN
            CALL D02NTF(NEQ,NEQMAX,JCEVAL,ML,MU,LENWJ,LENJP,WD(IPOSRW),
     *                  IFAIL1)
         ELSE IF (MATZ.EQ.'S') THEN
            CALL D02NUF(NEQ,NEQMAX,JCEVAL,LENWJ,IA,NIA,JA,NJA,IW(JACPOS)
     *                  ,LENJP,SENS,U,ETA,LBLOCK,ISPLIT,WD(IPOSRW),
     *                  IFAIL1)
         END IF
C
C  ...  Setup the ODE integrator  ...
C
C     OCT 94 -- PASS DEFAULT SNORM TO SETUP ROUTINES AND THEN ASSIGN WD(
C     TO CORRECT VALUE AFTERWARDS (TO ALLOW FOR OPTION OF L1 NORM FOR UP
C     SCHEME ROUTINES ONLY).
C
         TSNORM = 'D'
         IF (OPT(1).EQ.1.D0) THEN
            CALL D02NVF(NEQMAX,NY2DIM,MAXORD,METOD1,PETZLD,CONST,TCRIT,
     *                  HMIN,HMAX,H0,MAXSTP,MXHNIL,TSNORM,WD(IPOSRW),
     *                  IFAIL1)
C
         ELSE IF (OPT(1).EQ.2.D0) THEN
            CALL D02MWY(NEQMAX,NY2DIM,METOD1,METOD2,THETA,CONST,TCRIT,
     *                  HMIN,HMAX,H0,MAXSTP,MXHNIL,TSNORM,WD(IPOSRW),
     *                  IFAIL1)
         END IF
C
         TSNORM = SNORM
C        ENSURE UPPER CASE
         CALL E04UDU(TSNORM)
         IF (TSNORM.EQ.'1') THEN
            WD(IPOSRW+48) = 3.0D0
         ELSE IF (TSNORM.EQ.'2') THEN
            WD(IPOSRW+48) = 1.0D0
         END IF
C
      ELSE IF (IND.NE.1) THEN
         GO TO 280
      END IF
C
C    ... The RWORK array:
C
      IPOSRW = 1
C
C    ... The YDOTI array:
C
      IPOSYD = IPOSRW + LENRW
      LENYD = NEQMAX
C
C    ... The YSAVE array:
C
      IPOSYS = IPOSYD + LENYD
      LENYS = NY2DIM*NEQMAX
C
      ISPLIT = 0
      LBLOCK = .FALSE.
      INFPOS = 1
      JACPOS = 24
C
      LEWT = 51
      LACOR = LEWT + NEQMAX
      LSAVR = LACOR + NEQMAX
      IREVCM = 0
      LIWUSD = LENJP
      LRWUSD = LENWJ
      IIFAIL = 1
C
      IF (ITRACE.LE.0) THEN
         JTRACE = -1
      ELSE
         JTRACE = ITRACE
      END IF
   40 CONTINUE
      LTRACE = JTRACE
      CALL D02NNF(NEQ,NEQMAX,T,TOUT,Y,WD(IPOSYD),WD(IPOSRW),RTOL,ATOL,
     *            ITOL,IW(INFPOS),WD(IPOSYS),NY2DIM,WKJAC,LENWJ,
     *            IW(JACPOS),LENJP,IMON,INLN,IRES,IREVCM,LDERIV,ITASK,
     *            JTRACE,IIFAIL)
      LTRACE = ITRACE
C
      GO TO (80,60,80,100,120,80,100,140,160,
     *       180,80) IREVCM
C
      IW(5) = IW(6)
      GO TO 200
C
   60 CONTINUE
      CALL RESID(PDEFCN,BNDARY,DUMFCN,DUMBND,NEQ,WD(IPOSRW+TN-1),Y,
     *           WD(IPOSRW+LACOR-1),WD(IPOSRW+LSAVR-1),IRES,WKRES,
     *           NWKRES,PDEFN,FLXPFF,FLXPLF,BNDY,ODEFN)
C
      IF (IIFLAG.EQ.2) THEN
         GO TO 340
      ELSE
         GO TO 40
      END IF
C
   80 CONTINUE
C
      CALL RESID(PDEFCN,BNDARY,DUMFCN,DUMBND,NEQ,WD(IPOSRW+TN-1),Y,
     *           WD(IPOSYD),WD(IPOSRW+LSAVR-1),IRES,WKRES,NWKRES,PDEFN,
     *           FLXPFF,FLXPLF,BNDY,ODEFN)
C
      IF (IIFLAG.EQ.2) THEN
         GO TO 340
      ELSE
         GO TO 40
      END IF
C
  100 CONTINUE
      CALL RESID(PDEFCN,BNDARY,DUMFCN,DUMBND,NEQ,WD(IPOSRW+TN-1),Y,
     *           WD(IPOSYD),WD(IPOSRW+LACOR-1),IRES,WKRES,NWKRES,PDEFN,
     *           FLXPFF,FLXPLF,BNDY,ODEFN)
C
      IF (IIFLAG.EQ.2) THEN
         GO TO 340
      ELSE
         GO TO 40
      END IF
C
  120 CONTINUE
      CALL RESID(PDEFCN,BNDARY,DUMFCN,DUMBND,NEQ,WD(IPOSRW+TN-1),Y,
     *           WD(IPOSRW+LSAVR-1),WD(IPOSYD),IRES,WKRES,NWKRES,PDEFN,
     *           FLXPFF,FLXPLF,BNDY,ODEFN)
C
      IF (IIFLAG.EQ.2) THEN
         GO TO 340
      ELSE
         GO TO 40
      END IF
C
  140 CONTINUE
      IF (MATZ.EQ.'F') THEN
         CALL JACFUL(NEQ,WD(IPOSRW+TN-1),Y,WD(IPOSYD),WD(IPOSRW+H-1),
     *               WD(IPOSRW+EL0-1),WKJAC)
         GO TO 40
      ELSE IF (MATZ.EQ.'B') THEN
         CALL JACBND(NEQ,WD(IPOSRW+TN-1),Y,WD(IPOSYD),WD(IPOSRW+H-1),
     *               WD(IPOSRW+EL0-1),INT(WD(IPOSRW+ML-1)),
     *               INT(WD(IPOSRW+MU-1)),WKJAC)
         GO TO 40
C
      ELSE IF (MATZ.EQ.'S') THEN
         CALL D02NRF(J,IPLACE,IW(INFPOS))
C
         IF (IPLACE.EQ.1) THEN
            CALL JACSPS(NEQ,WD(IPOSRW+TN-1),Y,WD(IPOSYD),WD(IPOSRW+H-1),
     *                  WD(IPOSRW+EL0-1),J,WD(IPOSRW+LSAVR-1))
         ELSE
            ERRMSG =
     *' The routine attempted to call a routine to form the
     * Jacobian when a dummy was provided. '
            CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
            CALL JACSPS(NEQ,WD(IPOSRW+TN-1),Y,WD(IPOSYD),WD(IPOSRW+H-1),
     *                  WD(IPOSRW+EL0-1),J,WD(IPOSRW+LACOR-1))
         END IF
         GO TO 40
      END IF
C
  160 CONTINUE
      CALL MONITR(NEQ,WD(IPOSRW+TN-1),WD(IPOSRW+HU-1),WD(IPOSRW+H-1),Y,
     *            WD(IPOSYD),WD(IPOSYS),NEQMAX,WD(IPOSRW+LSAVR-1),
     *            WD(IPOSRW+LEWT-1),WKRES,NWKRES,WKMON,NWKMON,IMON,INLN,
     *            WD(IPOSRW+HMIN1-1),WD(IPOSRW+HMAX1-1),MONITF,IRES,
     *            IXFIX,NIXFIX)
C
C VP NOV94
      IF (IMON.EQ.-2) THEN
C        Error in remeshing routines
         IFAIL1 = 17
         GO TO 380
      END IF
C
  180 CONTINUE
      GO TO 40
  200 CONTINUE
C
C
      IF (IIFLAG.EQ.1) THEN
         ERRMSG =
     *' The diffusive term D was found to illegally depend on time
     * derivatives.'
         CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
         IFAIL1 = 14
         GO TO 380
      END IF
C
      IF (IIFAIL.NE.0) THEN
         IF (IIFAIL.EQ.1) THEN
            IFAIL1 = 9
            GO TO 380
         ELSE IF (IIFAIL.EQ.2) THEN
            IFAIL1 = 12
            GO TO 380
         ELSE IF (IIFAIL.EQ.3) THEN
            IFAIL1 = 2
            GO TO 380
         ELSE IF (IIFAIL.EQ.4 .OR. IIFAIL.EQ.5) THEN
            IFAIL1 = 3
            GO TO 380
         ELSE IF (IIFAIL.EQ.6) THEN
            IFAIL1 = 13
            GO TO 380
         ELSE IF (IIFAIL.EQ.7 .OR. IIFAIL.EQ.8) THEN
            IF (MATZ.EQ.'S') THEN
               ICALL = 1
               CALL D02NXF(ICALL,LIWREQ,LIWUSD,LRWREQ,LRWUSD,NLU,NNZ,
     *                     NGP,ISPLIT,IGROW,LBLOCK,NBLOCK,IW(INFPOS))
C
               IF ((LIWUSD.LT.LIWREQ) .OR. (LRWUSD.LT.LRWREQ)) THEN
C VP changed next few lines 15/9/94
                  IF ( .NOT. REMESH) THEN
                     LIWREQ = LIWREQ + (NIW-LIWUSD)
                     LIWUSD = NIW
                  ELSE
                     LIWREQ = LIWREQ + (NIW-LIWUSD) + 1 + NIXFIX
                     LIWUSD = NIW + 1 + NIXFIX
                  END IF
                  LRWREQ = LRWREQ + (NW-LRWUSD)
                  LRWUSD = NW
                  IFAIL1 = 15
                  GO TO 360
               ELSE
                  IFAIL1 = 4
               END IF
            ELSE
               IFAIL1 = 4
               GO TO 380
            END IF
         ELSE IF (IIFAIL.EQ.9) THEN
            IFAIL1 = 5
            GO TO 380
         ELSE IF (IIFAIL.EQ.10) THEN
            IFAIL1 = 11
            GO TO 380
         ELSE IF (IIFAIL.EQ.11) THEN
            IFAIL1 = 6
            GO TO 380
         ELSE IF (IIFAIL.EQ.13) THEN
            IFAIL1 = 10
            GO TO 380
         ELSE IF (IIFAIL.EQ.14) THEN
            IFAIL1 = 7
            GO TO 380
         END IF
      END IF
      GO TO 380
C
  220 CONTINUE
      ERRMSG =
     *' Routine entered with ALGOPT(3) (=R1) when ALGOPT(3) =
     *  1.0 for Newton method and 2.0 for f/iter.
     *  check array sizes. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,OPT(3),0.D0)
      IFAIL1 = 1
      GO TO 380
C
  240 CONTINUE
      ERRMSG =
     *' Routine entered with ALGOPT(4) (=R1) when ALGOPT(4) =
     *  2.0 to use the Petzold error test and 1.0 for not
     *  using it. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,OPT(4),0.D0)
      IFAIL1 = 1
      GO TO 380
C
  260 CONTINUE
      ERRMSG =
     *' Routine entered with ALGOPT(7) (=R1) when ALGOPT(7) = 1.0
     *  for switching between iteration methods and 2.0
     *  for not switching. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,OPT(7),0.D0)
      IFAIL1 = 1
      GO TO 380
C
  280 CONTINUE
      ERRMSG =
     *' Routine entered with IND(=I1), when IND = 0
     *  for a first call or restart and IND = 1, for
     *  continuing integration. '
      CALL D02NNQ(ERRMSG,1,1,IND,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 380
C
  300 CONTINUE
      ERRMSG =
     *' The real workspace is not big enough
     *  it is (=I1), it needs to be at least (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NWD,ITOTW,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 380
C
  320 CONTINUE
      ERRMSG =
     *' The Jacobian workspace is not big enough it
     *  is (=I1), it needs to be at least (=I2). '
      CALL D02NNQ(ERRMSG,1,2,LENWJ,I,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 380
C
  340 CONTINUE
      ERRMSG =
     *' One of the user supplied routines sets
     *  IRES (=I1) repeatedly. Integration can not start.'
      CALL D02NNQ(ERRMSG,1,1,IRES,0,0,0.0D0,0.0D0)
      IFAIL1 = 8
      GO TO 380
C
  360 CONTINUE
      ERRMSG =
     *' When using the sparsity option, either the required size of
     * integer workspace IW is (=I1) where you have used (=I2),'
      CALL D02NNQ(ERRMSG,1,2,LIWREQ,LIWUSD,0,0.0D0,0.0D0)
C
      ERRMSG =
     *' or the required size of real workspace W is (=I1) where you
     *  have used (=I2).'
      CALL D02NNQ(ERRMSG,1,2,LRWREQ,LRWUSD,0,0.0D0,0.0D0)
C
  380 CONTINUE
      RETURN
      END
