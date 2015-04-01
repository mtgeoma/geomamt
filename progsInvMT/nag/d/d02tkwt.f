      SUBROUTINE D02TKW(XI,XIOLD,Z,DMZ,RHS,DELZ,DELDMZ,DQZ,DQDMZ,G,W,V,
     *                  DGZ,DF,F,VALSTR,SLOPE,SCALE,DSCALE,ACCUM,UHIGH,
     *                  IPVTG,INTEGS,IPVTW,NFXPNT,FIXPNT,IFLAG,FFUN,
     *                  FJAC,GAFUN,GAJAC,GBFUN,GBJAC,GUESS,M,MAXORD,NEQ,
     *                  MSTAR,MSHINF,ZVAL,ZVAL1,RHO,COEF,KCOL,B,ACOL,
     *                  ASAVE,ALEFT,ARIGHT,NLBC,KD,N,NOLD,NMAX,NZ,NDMZ,
     *                  NONLIN,ICARE,IGUESS,NTOL,TOLIN,ROOT,WGTMSH,
     *                  WGTERR,JTOL,LTOL,ERR,ERREST,DGR,LDGR,ERMX,IERMX,
     *                  IJERMX,SX,SZ,SDMZ,SVMESH,IOUT,IPRINT)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose
C
C      This subroutine is the actual driver.  The nonlinear iteration
C      strategy is controlled here ( see [4] ). Upon convergence, ERRCHK
C      is called to test for satisfaction of the requested tolerances.
C
C   Arguments
C
C      XI     - current mesh
C      XIOLD  - previous mesh
C      Z      - solution values on each subinterval and boundary
C               conditions
C      DMZ    - high order solution derivatives at each collocation
C               point of each interval for each equation
C      RHS    - actual right hand side
C      DELZ   - solution of linear system (hence, correction to Z)
C      DELDMZ - uncondensed solution (hence, correction ot DMZ)
C      DQZ    - correction to Z
C      DQDMZ  - correction ot DMZ
C      G      - condensed matrix and its decomposition
C      W      - submatrices used in condensation for each interval
C      V      - submatrices used in condensation for each interval
C      DGZ
C      DF     - array for Jacobian of Ordinary differential equations
C      F      - array to hold ODE evaluations
C      VALSTR - array used to store solution approximations for use
C               in extrapolation error test
C      SLOPE  - array used in remeshing
C      SCALE  - scalinmg factors for Z
C      DSCALE - scaling factors for DMZ
C      ACCUM  - array used in remeshing
C      UHIGH  - array used in remeshing
C      IPVTG  - pivot index for decompostion of G
C      INTEGS - array to describe structure of linear system in G
C      IPVTW  - pivot index for each submatrix W
C      NFXPNT - number of fixed points in the mesh (<=0)
C      FIXPNT - array containing fixed mesh points
C      IFLAG  - indicates computation finished with
C               = -2  Newton iteration failed to converge
C               = -1  NMAX Too small to continue
C               =  0  singular matrix (hence almost certain that
C                     Jacobian formulated incorrectly)
C               =  1  successful solution
C      FFUN   - procedure to evaluate ODEs
C      FJAC  - procedure to evaluate Jacobian of ODEs
C      GSUB   - procedure to evaluate boundary conditions
C      DGSUB  - procedure to evaluate Jacobian of boundary conditions
C      GUESS  - procedure to evaluate initial solution (if supplied)
C      M      - orders of ODEs
C      MAXORD - maximal order of ODEs
C      NEQ    - number of ODEs
C      MSTAR  - number of solution components (sum(M(i),i=1,NEQ))
C      MSHINF - mesh control parameters
C      ZVAL   - array to hold approximation of solution (length MSTAR)
C      ZVAL1  - array to hold approximation of solution (NEQ x MAXORD)
C      RHO    - collocation points
C      COEF   - mesh dependent monomial coefficients
C      KCOL   - number of collocation points
C      B      - rk-basis coefficients for each interval
C      ACOL   - rk-basis coefficients for collocation points
C      ASAVE
C      ALEFT  - left hand end of mesh
C      ARIGHT - right hand end of mesh
C      NLBC   - number of boundary conditions at ALEFT
C      KD     - dimension of workspaces (NEQ*KCOL)
C      N      - number of points in current mesh
C      NOLD   - number of points in previous mesh
C      NMAX   - maximum number of meshpoints
C      NZ     - number of elements of Z, (n+1)*mstar
C      NDMZ   - number od elements of DMZ, n*kcol*neq
C      NONLIN - nonlinear=0, linear=1
C      ICARE  - =-1  no convergence occurred (used for regular problems)
C               = 0  a regular problem
C               = 1  a sensitive problem
C               = 2  used for continuation (see description of ipar(10)
C                    in colnew).
C      IGUESS - determines whether or not GUESS is called
C      NTOL   - number of tolereances
C      TOLIN  - values of tolerances
C      ROOT   - exponent for use in remeshing
C      WGTMSH - weights for use in remeshing
C      WGTERR - weights for use in extrapolation error test
C      JTOL   - index of tolerenaces with respect to MSTAR
C      LTOL   - index of tolereances with respect to NEQ
C      ERMX   - the maximum value of the error measure
C                 (<tolin(ijermx) for success)
C      IERMX  - index of subinterval where ermx first occurs
C      IJERMX - component of solution corresponding to ERMX
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine CONTRL)
C
C   Local variables of interest
C      CHECK  - maximum tolerance value, used as part of criteria for
C               checking for nonlinear iteration convergence
C      RELAX  - the relaxation factor for damped Newton iteration
C      RELMIN - minimum allowable value for RELAX  (otherwise the
C               Jacobian is considered singular).
C      RLXOLD - previous RELAX
C      RSTART - initial value for RELAX when problem is sensitive
C      IFRZ   - number of fixed Jacobian iterations
C      LMTFRZ - maximum value for IFRZ before performing a reinversion
C      ITER   - number of iterations (counted only when Jacobian
C               reinversions are performed).
C      IPRED  - = 0  if RELAX is determined by a correction
C               = 1  if RELAX is determined by a prediction
C      IFREEZ - = 0  if the Jacobian is to be updated
C               = 1  if the Jacobian is currently fixed (frozen)
C      ICONV  - = 0  if no previous convergence has been obtained
C               = 1  if convergence on a previous mesh has been obtained
C      RNORM  - norm of RHS (right hand side) for current iteration
C      RNOLD  - norm of RHS for previous iteration
C      ANSCL  - scaled norm of newton correction
C      ANFIX  - scaled norm of newton correction at next step
C      ANORM  - scaled norm of a correction obtained with jacobian fixed
C      IMESH  - a control variable for subroutines NEWMSH and ERRCHK
C               = 1  the current mesh resulted from mesh selection
C                    or is the initial mesh.
C               = 2  the current mesh resulted from doubling the
C                    previous mesh
C
C **********************************************************************
C
C  Constants for control of nonlinear iteration
C
C
C
C     .. Parameters ..
      INTEGER           FLAG, NUMBER, MLIMIT, ALTER
      PARAMETER         (FLAG=1,NUMBER=2,MLIMIT=3,ALTER=4)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALEFT, ARIGHT, ERMX
      INTEGER           ICARE, IERMX, IFLAG, IGUESS, IJERMX, IOUT,
     *                  IPRINT, KCOL, KD, LDGR, MAXORD, MSTAR, N, NDMZ,
     *                  NEQ, NFXPNT, NLBC, NMAX, NOLD, NONLIN, NTOL, NZ
      LOGICAL           SVMESH
C     .. Array Arguments ..
      DOUBLE PRECISION  ACCUM(NMAX+1), ACOL(28,7), ASAVE(28,4), B(28),
     *                  COEF(KCOL,KCOL), DELDMZ(NMAX*KCOL*NEQ),
     *                  DELZ(MSTAR*(NMAX+1)), DF(NEQ,NEQ,MAXORD),
     *                  DGR(LDGR,NEQ,MAXORD), DGZ(MSTAR),
     *                  DMZ(NMAX*KCOL*NEQ), DQDMZ(NMAX*KCOL*NEQ),
     *                  DQZ(MSTAR*(NMAX+1)), DSCALE(NEQ*KCOL*NMAX),
     *                  ERR(MSTAR), ERREST(MSTAR), F(NEQ), FIXPNT(*),
     *                  G((NMAX*2+1)*MSTAR*MSTAR), RHO(KCOL),
     *                  RHS(NMAX*KCOL*NEQ+MSTAR), ROOT(NTOL),
     *                  SCALE(MSTAR*(NMAX+1)), SDMZ(NMAX*KCOL*NEQ),
     *                  SLOPE(NMAX+1), SX(NMAX+1), SZ(MSTAR*(NMAX+1)),
     *                  TOLIN(NTOL), UHIGH(NEQ,2), V(KD,MSTAR,NMAX),
     *                  VALSTR(MSTAR,4*NMAX), W(KD,KD,NMAX),
     *                  WGTERR(MSTAR), WGTMSH(NTOL), XI(NMAX+1),
     *                  XIOLD(NMAX+1), Z(MSTAR*(NMAX+1)), ZVAL(MSTAR),
     *                  ZVAL1(NEQ,MAXORD)
      INTEGER           INTEGS(3,NMAX+2), IPVTG(MSTAR*(NMAX+1)),
     *                  IPVTW(KD,NMAX), JTOL(NTOL), LTOL(NTOL), M(NEQ),
     *                  MSHINF(4)
C     .. Subroutine Arguments ..
      EXTERNAL          FFUN, FJAC, GAFUN, GAJAC, GBFUN, GBJAC, GUESS
C     .. Local Scalars ..
      DOUBLE PRECISION  ANDIF, ANFIX, ANORM, ANSCL, ARG, CHECK, FACT,
     *                  FACTOR, PRECIS, RELAX, RELMIN, RLXOLD, RNOLD,
     *                  RNORM, RSTART
      INTEGER           I, ICONV, ICOR, IFIN, IFREEZ, IFRZ, IMESH, INZ,
     *                  IPRED, IT, ITER, IZ, J, K, LIMIT, LMTFRZ, MSING,
     *                  NOCONV, NP1, SN
      LOGICAL           SETRES, SVISIT, USER
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          D02TKJ, D02TKL, D02TKU, D02TKV, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SQRT
C     .. Executable Statements ..
C      save  svisit, sn
C
      USER = IGUESS .EQ. 1
      SVMESH = .FALSE.
      SN = 0
C
      PRECIS = 50.0D0*X02AJF()
      RELMIN = 1.D-3
      RSTART = 1.D-2
      LMTFRZ = 4
      LIMIT = 40
C
C Compute the maximum tolerance
C
      CHECK = 0.D0
      DO 20 I = 1, NTOL
         CHECK = MAX(TOLIN(I),CHECK)
   20 CONTINUE
      IMESH = 1
      ICONV = 0
      IF (NONLIN.EQ.0) ICONV = 1
      ICOR = 0
      NOCONV = 0
      MSING = 0
C
C The main iteration begins here.
C Loop here until error tolerances are satisfied or
C the code fails (due to a singular matrix or storage limitations)
C
   40 CONTINUE
C
C Initialization for a new mesh
C
      ITER = 0
      SETRES = NONLIN .NE. 0 .AND. ITER .EQ. 0
      IF (NONLIN.GT.0) GO TO 100
C
C The linear case.
C Set up and solve equations
C
      CALL D02TKV(MSING,XI,XIOLD,DUMMY,DUMMY,Z,DMZ,G,W,V,RHS,DUMMY,
     *            INTEGS,IPVTG,IPVTW,RNORM,0,FFUN,FJAC,GAFUN,GAJAC,
     *            GBFUN,GBJAC,GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,DGZ,DF,
     *            F,KD,MSTAR,NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,NOLD,
     *            NMAX,ZVAL,ZVAL1)
C
C Check for a singular matrix
C
      IF (MSING.EQ.0) GO TO 820
   60 CONTINUE
      IF (MSING.LT.0) GO TO 80
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a)')
     *     'DEB A local elimination matrix is singular'
         CALL X04BAF(IOUT,REC)
      END IF
      GO TO 1000
   80 CONTINUE
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a)') 'DEB The global BVP-matrix is singular'
         CALL X04BAF(IOUT,REC)
      END IF
      IFLAG = 0
      RETURN
C
C Iteration loop for nonlinear case
C Define the initial relaxation parameter (= relax)
C
  100 CONTINUE
      RELAX = 1.D0
C
C Check for previous convergence and problem sensitivity
C
      IF (ICARE.EQ.1 .OR. ICARE.EQ.(-1)) RELAX = RSTART
      IF (ICONV.EQ.0) GO TO 360
C
C Convergence on a previous mesh has been obtained. Thus
C we have a very good initial approximation for the Newton
C process.    Proceed with one full newton and then iterate
C with a fixed Jacobian.
C
      IFREEZ = 0
C
C Evaluate right hand side and its norm  and
C find the first Newton correction
C
      CALL D02TKV(MSING,XI,XIOLD,Z,DMZ,DELZ,DELDMZ,G,W,V,RHS,DQDMZ,
     *            INTEGS,IPVTG,IPVTW,RNOLD,1,FFUN,FJAC,GAFUN,GAJAC,
     *            GBFUN,GBJAC,GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,DGZ,DF,
     *            F,KD,MSTAR,NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,NOLD,
     *            NMAX,ZVAL,ZVAL1)
C
      IF (IPRINT.LT.0) THEN
         WRITE (REC,FMT='(a)') 'DEB Fixed Jacobian iteration'
         CALL X04BAF(IOUT,REC)
         WRITE (REC,FMT='(a,i2,a,e10.2)') 'DEB Iteration = ', ITER,
     *     '  Norm (rhs) = ', RNOLD
         CALL X04BAF(IOUT,REC)
      END IF
      GO TO 140
C
C Solve for the next iterate.
C The value of IFREEZ determines whether this is a full
C Newton step (=0) or a fixed jacobian iteration (=1).
C
  120 CONTINUE
      IF (IPRINT.LT.0) THEN
         WRITE (REC,FMT='(a,i2,a,e10.2)') 'DEB Iteration = ', ITER,
     *     '  Norm (rhs) = ', RNOLD
         CALL X04BAF(IOUT,REC)
      END IF
      RNOLD = RNORM
      SETRES = NONLIN .NE. 0 .AND. ITER .EQ. 0
      CALL D02TKV(MSING,XI,XIOLD,Z,DMZ,DELZ,DELDMZ,G,W,V,RHS,DUMMY,
     *            INTEGS,IPVTG,IPVTW,RNORM,3+IFREEZ,FFUN,FJAC,GAFUN,
     *            GAJAC,GBFUN,GBJAC,GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,
     *            DGZ,DF,F,KD,MSTAR,NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,
     *            NOLD,NMAX,ZVAL,ZVAL1)
C
C Check for a singular matrix
C
  140 CONTINUE
      IF (MSING.NE.0) GO TO 60
      IF (IFREEZ.EQ.1) GO TO 160
C
C A full newton step
C
      ITER = ITER + 1
      IFRZ = 0
  160 CONTINUE
C
C Update Z and DMZ, compute new RHS and its norm
C
      DO 180 I = 1, NZ
         Z(I) = Z(I) + DELZ(I)
  180 CONTINUE
      DO 200 I = 1, NDMZ
         DMZ(I) = DMZ(I) + DELDMZ(I)
  200 CONTINUE
      SETRES = NONLIN .NE. 0 .AND. ITER .EQ. 0
      CALL D02TKV(MSING,XI,XIOLD,Z,DMZ,DELZ,DELDMZ,G,W,V,RHS,DUMMY,
     *            INTEGS,IPVTG,IPVTW,RNORM,2,FFUN,FJAC,GAFUN,GAJAC,
     *            GBFUN,GBJAC,GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,DGZ,DF,
     *            F,KD,MSTAR,NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,NOLD,
     *            NMAX,ZVAL,ZVAL1)
C
C Check monotonicity. If the norm of  RHS  gets smaller,
C proceed with a fixed Jacobian; else proceed cautiously,
C as if convergence has not been obtained before (ICONV=0).
C
      IF (RNORM.LT.PRECIS) GO TO 800
      IF (RNORM.GT.RNOLD) GO TO 280
      IF (IFREEZ.EQ.1) GO TO 220
      IFREEZ = 1
      GO TO 120
C
C Verify that the linear convergence with fixed Jacobian
C is fast enough.
C
  220 CONTINUE
      IFRZ = IFRZ + 1
      IF (IFRZ.GE.LMTFRZ) IFREEZ = 0
      IF (RNOLD.LT.4.D0*RNORM) IFREEZ = 0
C
C...       check convergence (iconv = 1).
C
      DO 260 IT = 1, NTOL
         INZ = LTOL(IT)
         DO 240 IZ = INZ, NZ, MSTAR
            IF (ABS(DELZ(IZ)).GT.TOLIN(IT)*(ABS(Z(IZ))+1.D0)) GO TO 120
  240    CONTINUE
  260 CONTINUE
C
C...       convergence obtained
C
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a,i2,a)') 'DEB Convergence after ', ITER,
     *     ' iterations'
         CALL X04BAF(IOUT,REC)
      END IF
      GO TO 820
C
C...      convergence of fixed jacobian iteration failed.
C
  280 CONTINUE
      IF (IPRINT.LT.0) THEN
         WRITE (REC,FMT='(a,i2,a,e10.2)') 'DEB Iteration = ', ITER,
     *     '  Norm (rhs) = ', RNOLD
         CALL X04BAF(IOUT,REC)
         WRITE (REC,FMT='(a)') 'DEB Switch to damped Newton iteration'
         CALL X04BAF(IOUT,REC)
      END IF
      ICONV = 0
      RELAX = RSTART
      DO 300 I = 1, NZ
         Z(I) = Z(I) - DELZ(I)
  300 CONTINUE
      DO 320 I = 1, NDMZ
         DMZ(I) = DMZ(I) - DELDMZ(I)
  320 CONTINUE
C
C Update old mesh
C
      NP1 = N + 1
      DO 340 I = 1, NP1
         XIOLD(I) = XI(I)
  340 CONTINUE
      NOLD = N
C
      ITER = 0
C
C No previous convergence has been obtained. Proceed
C with the damped Newton method.
C Evaluate RHS and find the first Newton correction.
C
  360 CONTINUE
C
C Switch back to user supplied GUESS if no convergence
C
      IF (USER .AND. .NOT. SVMESH) IGUESS = 1
      SETRES = NONLIN .NE. 0 .AND. ITER .EQ. 0
      IF (IPRINT.LT.0) THEN
         WRITE (REC,FMT='(a)') 'DEB Full damped Newton iteration'
         CALL X04BAF(IOUT,REC)
      END IF
      CALL D02TKV(MSING,XI,XIOLD,Z,DMZ,DELZ,DELDMZ,G,W,V,RHS,DQDMZ,
     *            INTEGS,IPVTG,IPVTW,RNOLD,1,FFUN,FJAC,GAFUN,GAJAC,
     *            GBFUN,GBJAC,GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,DGZ,DF,
     *            F,KD,MSTAR,NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,NOLD,
     *            NMAX,ZVAL,ZVAL1)
C
C Check for a singular matrix
C
      IF (MSING.NE.0) GO TO 60
C
C Bookkeeping for first mesh
C
      IF (IGUESS.EQ.1) IGUESS = 0
C
C Find initial scaling
C
      CALL D02TKJ(N,MSTAR,KCOL,Z,XI,SCALE,DSCALE,M,NEQ,MAXORD)
      GO TO 480
C
C Main iteration loop
C
  380 CONTINUE
      RNOLD = RNORM
      IF (ITER.GE.LIMIT) GO TO 940
C
C Update scaling
C
      CALL D02TKJ(N,MSTAR,KCOL,Z,XI,SCALE,DSCALE,M,NEQ,MAXORD)
C
C Compute norm of newton correction with new scaling
C
      ANSCL = 0.D0
      DO 400 I = 1, NZ
         ANSCL = ANSCL + (DELZ(I)*SCALE(I))**2
  400 CONTINUE
      DO 420 I = 1, NDMZ
         ANSCL = ANSCL + (DELDMZ(I)*DSCALE(I))**2
  420 CONTINUE
      ANSCL = SQRT(ANSCL/DBLE(NZ+NDMZ))
C
C Find a newton direction
C
      CALL D02TKV(MSING,XI,XIOLD,Z,DMZ,DELZ,DELDMZ,G,W,V,RHS,DUMMY,
     *            INTEGS,IPVTG,IPVTW,RNORM,3,FFUN,FJAC,GAFUN,GAJAC,
     *            GBFUN,GBJAC,GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,DGZ,DF,
     *            F,KD,MSTAR,NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,NOLD,
     *            NMAX,ZVAL,ZVAL1)
C
C Check for a singular matrix
C
      IF (MSING.NE.0) GO TO 60
C
C Predict relaxation factor for Newton step.
C
      ANDIF = 0.D0
      DO 440 I = 1, NZ
         ANDIF = ANDIF + ((DQZ(I)-DELZ(I))*SCALE(I))**2
  440 CONTINUE
      DO 460 I = 1, NDMZ
         ANDIF = ANDIF + ((DQDMZ(I)-DELDMZ(I))*DSCALE(I))**2
  460 CONTINUE
      ANDIF = SQRT(ANDIF/DBLE(NZ+NDMZ)+PRECIS)
      RELAX = RELAX*ANSCL/ANDIF
      IF (RELAX.GT.1.D0) RELAX = 1.D0
  480 CONTINUE
      RLXOLD = RELAX
      IPRED = 1
      ITER = ITER + 1
C
C Determine a new Z and DMZ and find new RHS and its norm
C
      DO 500 I = 1, NZ
         Z(I) = Z(I) + RELAX*DELZ(I)
  500 CONTINUE
      DO 520 I = 1, NDMZ
         DMZ(I) = DMZ(I) + RELAX*DELDMZ(I)
  520 CONTINUE
  540 CONTINUE
      SETRES = NONLIN .NE. 0 .AND. ITER .EQ. 0
      CALL D02TKV(MSING,XI,XIOLD,Z,DMZ,DQZ,DQDMZ,G,W,V,RHS,DUMMY,INTEGS,
     *            IPVTG,IPVTW,RNORM,2,FFUN,FJAC,GAFUN,GAJAC,GBFUN,GBJAC,
     *            GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,DGZ,DF,F,KD,MSTAR,
     *            NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,NOLD,NMAX,ZVAL,
     *            ZVAL1)
C
C Compute a fixed Jacobian iterate (used to control RELAX)
C
      CALL D02TKV(MSING,XI,XIOLD,Z,DMZ,DQZ,DQDMZ,G,W,V,RHS,DUMMY,INTEGS,
     *            IPVTG,IPVTW,RNORM,4,FFUN,FJAC,GAFUN,GAJAC,GBFUN,GBJAC,
     *            GUESS,IGUESS,M,MAXORD,NEQ,DGR,LDGR,DGZ,DF,F,KD,MSTAR,
     *            NLBC,SETRES,KCOL,RHO,COEF,B,ACOL,N,NOLD,NMAX,ZVAL,
     *            ZVAL1)
C
C Find scaled norms of various terms used to correct RELAX
C
      ANORM = 0.D0
      ANFIX = 0.D0
      DO 560 I = 1, NZ
         ANORM = ANORM + (DELZ(I)*SCALE(I))**2
         ANFIX = ANFIX + (DQZ(I)*SCALE(I))**2
  560 CONTINUE
      DO 580 I = 1, NDMZ
         ANORM = ANORM + (DELDMZ(I)*DSCALE(I))**2
         ANFIX = ANFIX + (DQDMZ(I)*DSCALE(I))**2
  580 CONTINUE
      ANORM = SQRT(ANORM/DBLE(NZ+NDMZ))
      ANFIX = SQRT(ANFIX/DBLE(NZ+NDMZ))
C
      IF (IPRINT.LT.0) THEN
         IF (ICOR.EQ.1) THEN
            WRITE (REC,FMT='(a,e10.2)')
     *        'DEB Relaxation factor corrected to ', RELAX
         ELSE
            WRITE (REC,FMT='(a,i2,a,e10.2)') 'DEB Iteration = ', ITER,
     *        '  Relaxation factor = ', RELAX
         END IF
         CALL X04BAF(IOUT,REC)
         WRITE (REC,FMT='(a,e10.2,a,e10.2)')
     *     'DEB Norm of scaled rhs changes from ', ANORM, ' to ', ANFIX
         CALL X04BAF(IOUT,REC)
         WRITE (REC,FMT='(a,e10.2,a,e10.2)')
     *     'DEB Norm  of   rhs   changes  from  ', RNOLD, ' to ', RNOLD
         CALL X04BAF(IOUT,REC)
      END IF
C      IF (ICOR.EQ.1) GO TO 600
C      IF (IPRINT.LT.0) WRITE (IOUT,FMT=99995) ITER, RELAX, ANORM, ANFIX
C     *    RNOLD, RNORM
C      GO TO 620
C  600 CONTINUE
C      IF (IPRINT.LT.0) WRITE (IOUT,FMT=99992) RELAX, ANORM, ANFIX,
C     *    RNOLD, RNORM
C  620 CONTINUE
      ICOR = 0
C
C Check for monotonic decrease in DELZ and DELDMZ.
C
      IF (ANFIX.LT.PRECIS .OR. RNORM.LT.PRECIS) GO TO 800
      IF (ANFIX.GT.ANORM) GO TO 600
C
C We have a decrease.
C If DQZ and DQDMZ are small, check for convergence
C
      IF (ANFIX.LE.CHECK) GO TO 700
C
C Correct the predicted RELAX unless the corrected
C value is within 10 percent of the predicted one.
C
      IF (IPRED.NE.1) GO TO 380
  600 CONTINUE
      IF (ITER.GE.LIMIT) GO TO 940
C
C Correct the relaxation factor.
C
      IPRED = 0
      ARG = (ANFIX/ANORM-1.D0)/RELAX + 1.D0
      IF (ARG.LT.0.D0) GO TO 380
      IF (ARG.LE..25D0*RELAX+.125D0*RELAX**2) GO TO 620
      FACTOR = -1.D0 + SQRT(1.D0+8.D0*ARG)
      IF (ABS(FACTOR-1.D0).LT..1D0*FACTOR) GO TO 380
      IF (FACTOR.LT.0.5D0) FACTOR = 0.5D0
      RELAX = RELAX/FACTOR
      GO TO 640
  620 CONTINUE
      IF (RELAX.GE.0.9D0) GO TO 380
      RELAX = 1.D0
  640 CONTINUE
      ICOR = 1
      IF (RELAX.LT.RELMIN) GO TO 960
      FACT = RELAX - RLXOLD
      DO 660 I = 1, NZ
         Z(I) = Z(I) + FACT*DELZ(I)
  660 CONTINUE
      DO 680 I = 1, NDMZ
         DMZ(I) = DMZ(I) + FACT*DELDMZ(I)
  680 CONTINUE
      RLXOLD = RELAX
      GO TO 540
C
C Check convergence (ICONV = 0).
C
  700 CONTINUE
      DO 740 IT = 1, NTOL
         INZ = LTOL(IT)
         DO 720 IZ = INZ, NZ, MSTAR
            IF (ABS(DQZ(IZ)).GT.TOLIN(IT)*(ABS(Z(IZ))+1.D0)) GO TO 380
  720    CONTINUE
  740 CONTINUE
C
C Convergence obtained
C
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a,i2,a)') 'DEB Convergence after ', ITER,
     *     ' iterations'
         CALL X04BAF(IOUT,REC)
      END IF
C
C Since convergence obtained, update Z and DMZ with term
C from the fixed Jacobian iteration.
C
      DO 760 I = 1, NZ
         Z(I) = Z(I) + DQZ(I)
  760 CONTINUE
      DO 780 I = 1, NDMZ
         DMZ(I) = DMZ(I) + DQDMZ(I)
  780 CONTINUE
  800 CONTINUE
      IF ((ANFIX.LT.PRECIS .OR. RNORM.LT.PRECIS) .AND. IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a,i2,a)') 'DEB Convergence after ', ITER,
     *     ' iterations'
         CALL X04BAF(IOUT,REC)
      END IF
      ICONV = 1
      IF (ICARE.EQ.(-1)) ICARE = 0
C
C If full output has been requested, print values of the
C solution components Z at the meshpoints.
C
  820 CONTINUE
      IF (IPRINT.LT.0) THEN
         DO 860 J = 1, MSTAR
            WRITE (REC,FMT='(a,i2,a)') 'DEB Mesh values for Z(', J, '):'
            CALL X04BAF(IOUT,REC)
            DO 840 I = J, NZ, 5*MSTAR
               WRITE (REC,FMT='(a,5e15.7)') 'DEB ',
     *           (Z(K),K=I,MIN(I+4*MSTAR,NZ),MSTAR)
               CALL X04BAF(IOUT,REC)
  840       CONTINUE
  860    CONTINUE
      END IF
C
C Save this solution in case we need to go back to it
C
      SN = N
C      snz = nz
C      sndmz = ndmz
      SVISIT = .FALSE.
      SVMESH = .TRUE.
      DO 880 I = 1, N + 1
         SX(I) = XI(I)
  880 CONTINUE
      DO 900 I = 1, NZ
         SZ(I) = Z(I)
  900 CONTINUE
      DO 920 I = 1, NDMZ
         SDMZ(I) = DMZ(I)
  920 CONTINUE
C
C Check for error tolerance satisfaction
C
      IFIN = 1
      IF (IMESH.EQ.2) THEN
         CALL D02TKU(XI,N,Z,DMZ,VALSTR,IFIN,KCOL,NEQ,M,MAXORD,MSTAR,
     *               ASAVE,ERR,ERREST,WGTERR,TOLIN,LTOL,NTOL,ERMX,IERMX,
     *               IJERMX)
         IF (IPRINT.LT.1) THEN
            IF (IFIN.EQ.1) THEN
               WRITE (REC,FMT='(a)') 'DEB Error test passed'
               CALL X04BAF(IOUT,REC)
            ELSE
               WRITE (REC,FMT='(a)') 'DEB Error test failed:'
               CALL X04BAF(IOUT,REC)
               WRITE (REC,FMT='(a,i6,a,i3)')
     *           'DEB    worst on subinterval ', IERMX,
     *           ' for component ', IJERMX
               CALL X04BAF(IOUT,REC)
               WRITE (REC,FMT='(a,e10.3,a,e10.3)')
     *           'DEB    with measure ', ERMX, ' > TOL = ',
     *           TOLIN(IJERMX)
               CALL X04BAF(IOUT,REC)
            END IF
         END IF
C
C MSHINF(FLAG) indicates that solution values on the current mesh have b
C computed which may then be used on the next error test if this current
C mesh is extrapolated upon.
C
         MSHINF(FLAG) = 1
      END IF
C
      IF (IMESH.EQ.1 .OR. IFIN.EQ.0 .AND. ICARE.NE.2) GO TO 1000
      IFLAG = 1
      RETURN
C
C Diagnostics for failure of nonlinear iteration.
C
  940 CONTINUE
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a,i2,a)') 'DEB No convergence after ', ITER,
     *     ' iterations'
         CALL X04BAF(IOUT,REC)
      END IF
      GO TO 980
  960 CONTINUE
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a,e10.3,a,a,e10.3,a)')
     *     'DEB No convergence. Relaxation factor =', RELAX, ' is too ',
     *     'small (<', RELMIN, ')'
         CALL X04BAF(IOUT,REC)
      END IF
  980 CONTINUE
      IFLAG = -2
      NOCONV = NOCONV + 1
      IF (ICARE.EQ.2 .AND. NOCONV.GT.1) THEN
         N = SN
         RETURN
      END IF
      IF (ICARE.EQ.0) ICARE = -1
C
C...       update old mesh
C
 1000 CONTINUE
C
C Try something different
C
      IF (ICONV.EQ.0 .AND. SVMESH) THEN
         NOLD = SN
C We're retreating to the solution on a previous mesh
C We're going to halve because we had a convergence failure
C The mesh that we halve is that associated with  the
C previously converged solution - if that hasn't already been doubled.
         IF (IPRINT.LT.1) THEN
            WRITE (REC,FMT='(a,i6,a)')
     *        'DEB Reverting to a converged solution on the previous ',
     *        SN + 1, ' point mesh.'
            CALL X04BAF(IOUT,REC)
         END IF
         IF ( .NOT. SVISIT) THEN
            IF (IPRINT.LT.1) THEN
               WRITE (REC,FMT='(a)')
     *           'DEB Will try doubling the number of mesh points.'
               CALL X04BAF(IOUT,REC)
            END IF
            ICONV = 1
            N = NOLD
            DO 1020 I = 1, NOLD + 1
               XI(I) = SX(I)
 1020       CONTINUE
         ELSE
            IF (IPRINT.LT.1) THEN
               WRITE (REC,FMT='(a)')
     *          'DEB Will try doubling the number of mesh points again.'
               CALL X04BAF(IOUT,REC)
            END IF
         END IF
         DO 1040 I = 1, NOLD + 1
            XIOLD(I) = SX(I)
 1040    CONTINUE
         NZ = MSTAR*(NOLD+1)
         DO 1060 I = 1, NZ
            Z(I) = SZ(I)
 1060    CONTINUE
         NDMZ = NEQ*KCOL*NOLD
         DO 1080 I = 1, NDMZ
            DMZ(I) = SDMZ(I)
 1080    CONTINUE
         IF (N.GT.NMAX/2) THEN
            N = 2*N
            GO TO 1120
         END IF
         SVISIT = .TRUE.
C
      ELSE
C
         NP1 = N + 1
         DO 1100 I = 1, NP1
            XIOLD(I) = XI(I)
 1100    CONTINUE
         NOLD = N
C
      END IF
C
C Pick a new mesh
C Check safeguards for mesh refinement
C
      IMESH = 1
      IF (ICONV.EQ.0 .OR. MSHINF(NUMBER).GE.MSHINF(MLIMIT)
     *    .OR. MSHINF(ALTER).GE.MSHINF(MLIMIT)) IMESH = 2
      IF (MSHINF(ALTER).GE.MSHINF(MLIMIT) .AND. MSHINF(NUMBER)
     *    .LT.MSHINF(MLIMIT)) MSHINF(ALTER) = 1

      CALL D02TKL(IMESH,ALEFT,ARIGHT,XI,XIOLD,N,NOLD,NMAX,Z,DMZ,VALSTR,
     *            SLOPE,ACCUM,UHIGH,NFXPNT,FIXPNT,KCOL,NEQ,M,MAXORD,
     *            MSTAR,COEF,IGUESS,ASAVE,NTOL,JTOL,LTOL,WGTMSH,ROOT,
     *            MSHINF,IOUT,IPRINT)
      NZ = MSTAR*(N+1)
      NDMZ = NEQ*KCOL*N
C
C Exit if expected N is too large (but may try N = NMAX once)
C
 1120 CONTINUE
      IF (N.GT.NMAX) THEN
         N = N/2
         IFLAG = -1
         IF (IPRINT.LT.1) THEN
            WRITE (REC,FMT='(a,i6)')
     *        'DEB Would exceed maximum number of mesh points ',
     *        NMAX + 1
            CALL X04BAF(IOUT,REC)
            WRITE (REC,FMT='(a,i6)') 'DEB Would want to try at least ',
     *        2*N
            CALL X04BAF(IOUT,REC)
         END IF
         IF (ICONV.EQ.0) THEN
            IF (IPRINT.LT.1) THEN
               WRITE (REC,FMT='(a)')
     *           'DEB Haven''t been able to converge despite doubling'
               CALL X04BAF(IOUT,REC)
            END IF
            IFLAG = -2
         END IF
         IF (ICONV.EQ.1 .AND. IPRINT.LT.1) THEN
            WRITE (REC,FMT='(a,a)')
     *        'DEB Could try relaxing tolerances or increasing max ',
     *        'number of points'
            CALL X04BAF(IOUT,REC)
         END IF
         N = SN
         RETURN
      END IF
C
      IF (ICONV.EQ.0) IMESH = 1
      IF (ICARE.EQ.1) ICONV = 0
      GO TO 40
C
      END
