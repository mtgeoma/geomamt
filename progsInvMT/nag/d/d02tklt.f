      SUBROUTINE D02TKL(MODE,ALEFT,ARIGHT,XI,XIOLD,N,NOLD,NMAX,Z,DMZ,
     *                  VALSTR,SLOPE,ACCUM,UHIGH,NFXPNT,FIXPNT,KCOL,NEQ,
     *                  M,MAXORD,MSTAR,COEF,IGUESS,ASAVE,NTOL,JTOL,LTOL,
     *                  WGTMSH,ROOT,MSHINF,IOUT,IPRINT)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C**********************************************************************
C   Purpose:
C      Select a mesh on which a collocation solution is to be
C      determined
C      There are 5 possible modes of action:
C         MODE = 5,4,3 - deal mainly with definition of an initial
C                        mesh for the current boundary value problem
C              = 2,1   - deal with definition of a new mesh, either
C                        by simple mesh halving or by mesh selection
C      More specifically, for
C         MODE = 5  an initial (generally nonuniform) mesh is
C                   defined by the user and no mesh selection is to
C                   be performed
C              = 4  an initial (generally nonuniform) mesh is
C                   defined by the user
C              = 3  a simple uniform mesh (except possibly for some
C                   fixed points) is defined; n= no. of subintervals
C              = 1  the automatic mesh selection procedure is used
C              = 2  a simple mesh halving is performed
C   Arguments:
C      MODE   - Action to be takne (see above)
C      ALEFT  - Leftmost point of mesh
C      ARIGHT - rightmost point of mesh
C      XI     - the mesh selected
C      XIOLD  - the previous mesh (if any)
C      N      - number of mesh subintervals in XI
C      NOLD   - number of subintervals for former mesh XIOLD
C      NMAX   - maximum number of subintervals permitted
C      Z      - current solution vector z(u(x)) for each interval
C      DMZ    - array of m(j)-th derivatives of the current
C               solution at each collocation point for each interval
C      VALSTR - set of solution vectors z(u(x)) sampled at four points
C               within each interval (used in the error test)
C      SLOPE  - estimate of a high order derivative of solution
C               on each interval for use in remeshing; an approximate
C               quantity to be redistributed
C                             .                        (k+m(j))
C                     slope(i)=     max   (weight(l) *u      (xi(i)))
C                               1 <= l <= ntol         j
C                     where j=jtol(l)
C      ACCUM  - the integral of SLOPE across the range
C               ACCUM(i) = integral from ALEFT to XIOLD(i)
C      UHIGH  - array used in computing high order derivative approx
C      NFXPNT - number of fixed points within the mesh
C      FIXPNT - array of fixed points within the mesh
C      KCOL   - number of collocation points per interval
C      NEQ    - number of ODEs
C      M      - orderrs of ODEs
C      MAXORD - maximal order of ODEs
C      MSTAR  - sum(M(i),i=1,NEQ), length of solution vector for
C               each interval
C      COEF   - rk-basis coefficients for each interval
C      IGUESS - specifies nature of initial mesh for MODE=5 and 4
C               = 2 the subroutine sets XIi=XIOLD;  this is used e.g. if
C                   continuation is being performed, and a mesh for the
C                   old differential equation is being used
C               = 3 same as for =2, except XI uses every other point of
C                   XIOLD (so mesh xiold is mesh XI halved)
C               = 4 XI has been defined by the user, and an old mesh
C                   XIOLD is also available otherwise, XI has been
C                   defined by the user and we set XIOLD=XI here
C      ASAVE  - rk-basis coefficients for sampling points in each
c               interval used in the error test
C      NTOL   - number of solution components used in error measure
C               and control
C      JTOL   - indices of ODEs which  used in error measure and control
C      LTOL   - indices of components in Z used in erro measure and
C               control
C      WGTMSH - scaling factors used in error measure
C      ROOT   - exponents used in error measure
C      MSHINF - specifies info about the mesh selection strategy
C               (1=FLAG)
C                        = 1  the mesh is a halving of its former mesh
C                            (so an error estimate has been calculated)
C                        = 0  otherwise
C               (2=NUMBER)
C                        =  number of mesh selections which have
C                           actually been performed for the given N
C               (3=MLIMIT)
C                        = maximum number of mesh selections which are
C                          permitted for a given N before mesh halving
C               (4=ALTER)
C                    = number of consecutive times ( plus 1 ) the mesh
C                      selection has alternately halved and doubled N;
C                      if >= MSHINF(MLIMIT) then  CONTRL  requires
C                      that the current mesh be halved.
C      IOUT   - channel for any trace information
C      IPRINT - trace indicator
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine NEWMSH)
C**********************************************************************
C     .. Parameters ..
      INTEGER           FLAG, NUMBER, MLIMIT, ALTER
      PARAMETER         (FLAG=1,NUMBER=2,MLIMIT=3,ALTER=4)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALEFT, ARIGHT
      INTEGER           IGUESS, IOUT, IPRINT, KCOL, MAXORD, MODE, MSTAR,
     *                  N, NEQ, NFXPNT, NMAX, NOLD, NTOL
C     .. Array Arguments ..
      DOUBLE PRECISION  ACCUM(NMAX+1), ASAVE(28,4), COEF(KCOL,KCOL),
     *                  DMZ(NEQ,KCOL,NMAX), FIXPNT(*), ROOT(NTOL),
     *                  SLOPE(NMAX), UHIGH(NEQ,2), VALSTR(MSTAR,4*NMAX),
     *                  WGTMSH(NTOL), XI(NMAX+1), XIOLD(NMAX+1),
     *                  Z(MSTAR,NMAX)
      INTEGER           JTOL(NTOL), LTOL(NTOL), M(NEQ), MSHINF(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACCL, ACCR, AVRG, DEGEQU, DX, HD6, HIOLD,
     *                  ONEOVH, SLPHMX, TEMP, TSUM, X, XLEFT, XRIGHT
      INTEGER           I, IFLIP, ILEFT, IN, IOLD, IRIGHT, J, JJ, JZ,
     *                  KSTORE, L, LCARRY, LNEW, LOLD, N2, NACCUM,
     *                  NFXP1, NMAX2, NMIN, NMX, NOLDP1, NP1, NREGN
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           D02TKS
      EXTERNAL          X02AJF, D02TKS
C     .. External Subroutines ..
      EXTERNAL          D02TKK, D02TKT, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN
C     .. Executable Statements ..
      NFXP1 = NFXPNT + 1
      GO TO (260,160,100,40,20) MODE
C 
C  MODE = 5 - set MSHINF(MLIMIT)=1 so that no mesh selection is performe
C 
   20 CONTINUE
      MSHINF(MLIMIT) = 1
C 
C  MODE = 4 - the user-specified initial mesh is already in place.
C 
   40 CONTINUE
      IF (IGUESS.GE.2) THEN
C 
C  IGUESS = 2, 3 or 4.
C 
         NOLDP1 = NOLD + 1
         IF (IPRINT.LT.1) THEN
            WRITE (REC,FMT='(a,i6,a)') 'DEB Former mesh of ', NOLDP1,
     *        ' points:'
            CALL X04BAF(IOUT,REC)
            DO 60 I = 1, NOLDP1, 5
               WRITE (REC,FMT='(a,5e14.6)') 'DEB ',
     *           (XIOLD(J),J=I,MIN(I+4,NOLDP1))
               CALL X04BAF(IOUT,REC)
   60       CONTINUE
         END IF
         IF (IGUESS.EQ.3) THEN
C 
C  If IREAD ( ipar(8) ) .ge. 1 and IGUESS ( ipar(9) ) .eq. 3
C  then the first mesh is every second point of the
C  mesh in  XIOLD.
C 
            N = NOLD/2
            I = 0
            DO 80 J = 1, NOLD, 2
               I = I + 1
               XI(I) = XIOLD(J)
   80       CONTINUE
         END IF
      END IF
      NP1 = N + 1
      XI(1) = ALEFT
      XI(NP1) = ARIGHT
      GO TO 480
C 
C  MODE = 3 -  generate a (piecewise) uniform mesh. If there are
C  fixed points then ensure that the n being used is large enough.
C 
  100 CONTINUE
      IF (N.LT.NFXP1) N = NFXP1
      NP1 = N + 1
      XI(1) = ALEFT
      ILEFT = 1
      XLEFT = ALEFT
C 
C  Loop over the subregions between fixed points.
C 
      DO 140 J = 1, NFXP1
         XRIGHT = ARIGHT
         IRIGHT = NP1
         IF (J.LT.NFXP1) THEN
            XRIGHT = FIXPNT(J)
C 
C  Determine where the J-th fixed point should fall in the
C  new mesh - this is XI(IRIGHT) and the (J-1)st fixed
C  point is in XI(ILEFT).
C 
            NMIN = (XRIGHT-ALEFT)/(ARIGHT-ALEFT)*DBLE(N) + 1.5D0
            IF (NMIN.GT.N-NFXPNT+J) NMIN = N - NFXPNT + J
            IRIGHT = MAX(ILEFT+1,NMIN)
         END IF
         XI(IRIGHT) = XRIGHT
C 
C  Generate equally spaced points between the (J-1)st and the
C  J-th fixed points.
C 
         NREGN = IRIGHT - ILEFT - 1
         IF (NREGN.GT.0) THEN
            DX = (XRIGHT-XLEFT)/DBLE(NREGN+1)
            DO 120 I = 1, NREGN
               XI(ILEFT+I) = XLEFT + DBLE(I)*DX
  120       CONTINUE
         END IF
         ILEFT = IRIGHT
         XLEFT = XRIGHT
  140 CONTINUE
      GO TO 480
C 
C  MODE = 2 - halve the current mesh (i.e. double its size).
C 
  160 CONTINUE
      N2 = 2*N
C 
C  Check that N does not exceed storage limitations.
C 
      IF (N2.GT.NMAX) THEN
C 
C  If possible, try with N = NMAX. Redistribute first.
C 
         IF (MODE.NE.2) THEN
            N = NMAX/2
            GO TO 340
         ELSE
            IF (IPRINT.LT.1) THEN
               WRITE (REC,FMT='(a,i6,a)')
     *           'DEB Expected number of mesh points ', N2 - 1,
     *           ' too large'
               CALL X04BAF(IOUT,REC)
            END IF
            N = N2
            RETURN
         END IF
      END IF
C 
C  Calculate the old approximate solution values at points to be used in
C  ERRCHK for error estimates. If MSHINF(FLAG) = 1 an error estimate was
C  obtained for for the old approximation so half the needed values will
C  already be in VALSTR.
C 
      IF (MSHINF(FLAG).EQ.1) THEN
C 
C  Save in VALSTR the values of the old solution
C  at the relative positions 1/6 and 5/6 in each subinterval.
C 
         KSTORE = 1
         DO 180 I = 1, NOLD
            HD6 = (XIOLD(I+1)-XIOLD(I))/6.0D0
            X = XIOLD(I) + HD6
            CALL D02TKT(X,VALSTR(1,KSTORE),ASAVE(1,1),DUMMY,XIOLD(I),
     *                  XIOLD(I+1),Z(1,I),DMZ(1,1,I),KCOL,NEQ,MAXORD,M,
     *                  MSTAR,.FALSE.,.FALSE.,DUMMY)
            X = X + 4.0D0*HD6
            KSTORE = KSTORE + 3
            CALL D02TKT(X,VALSTR(1,KSTORE),ASAVE(1,4),DUMMY,XIOLD(I),
     *                  XIOLD(I+1),Z(1,I),DMZ(1,1,I),KCOL,NEQ,MAXORD,M,
     *                  MSTAR,.FALSE.,.FALSE.,DUMMY)
            KSTORE = KSTORE + 1
  180    CONTINUE
      ELSE
C 
C  Save in  VALSTR  the values of the old solution
C  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
C  each subinterval.
C 
         KSTORE = 1
         DO 220 I = 1, N
            X = XI(I)
            HD6 = (XI(I+1)-XI(I))/6.0D0
            DO 200 J = 1, 4
               X = X + HD6
               IF (J.EQ.3) X = X + HD6
               IOLD = D02TKS(X,XIOLD,NOLD)
               CALL D02TKT(X,VALSTR(1,KSTORE),ASAVE(1,J),DUMMY,
     *                     XIOLD(IOLD),XIOLD(IOLD+1),Z(1,IOLD),
     *                     DMZ(1,1,IOLD),KCOL,NEQ,MAXORD,M,MSTAR,
     *                     .FALSE.,.FALSE.,DUMMY)
               KSTORE = KSTORE + 1
  200       CONTINUE
  220    CONTINUE
      END IF
C 
      MSHINF(FLAG) = 0
      MSHINF(NUMBER) = 1
      MODE = 2
C 
C  Generate the halved mesh. Halve XI - note that XIOLD might contain
C  a converged mesh (not necessarily the previous XI!).
C 
      J = 2*N
      DO 240 I = N, 1, -1
         XI(J+1) = XI(I+1)
         XI(J) = (XI(I)+XI(I+1))/2.0D0
         J = J - 2
  240 CONTINUE
      N = N2
      GO TO 480
C 
C  MODE = 1  we do mesh selection if it is deemed worthwhile.
C 
  260 CONTINUE
      IF (NOLD.EQ.1) GO TO 160
      IF (NOLD.LE.2*NFXPNT) GO TO 160
C 
C  The first interval has to be treated separately from the
C  other intervals (generally the solution on the (I-1)st and ith
C  intervals will be used to approximate the needed derivative, but
C  here the 1st and second intervals are used.)
C 
      I = 1
      HIOLD = XIOLD(2) - XIOLD(1)
      CALL D02TKK(UHIGH(1,1),HIOLD,DMZ(1,1,1),NEQ,KCOL,COEF)
      HIOLD = XIOLD(3) - XIOLD(2)
      CALL D02TKK(UHIGH(1,2),HIOLD,DMZ(1,1,2),NEQ,KCOL,COEF)
      ACCUM(1) = 0.0D0
      SLOPE(1) = 0.0D0
      ONEOVH = 2.0D0/(XIOLD(3)-XIOLD(1))
      DO 280 J = 1, NTOL
         JJ = JTOL(J)
         JZ = LTOL(J)
         SLOPE(1) = MAX(SLOPE(1),(ABS(UHIGH(JJ,2)-UHIGH(JJ,1))*WGTMSH(J)
     *              *ONEOVH/(1.0D0+ABS(Z(JZ,1))))**ROOT(J))
  280 CONTINUE
      SLPHMX = SLOPE(1)*(XIOLD(2)-XIOLD(1))
      ACCUM(2) = SLPHMX
      IFLIP = 1
C 
C  Go through the remaining intervals generating SLOPE and  ACCUM.
C 
      DO 320 I = 2, NOLD
         HIOLD = XIOLD(I+1) - XIOLD(I)
         IF (IFLIP.EQ.-1) CALL D02TKK(UHIGH(1,1),HIOLD,DMZ(1,1,I),NEQ,
     *                                KCOL,COEF)
         IF (IFLIP.EQ.1) CALL D02TKK(UHIGH(1,2),HIOLD,DMZ(1,1,I),NEQ,
     *                               KCOL,COEF)
         ONEOVH = 2.0D0/(XIOLD(I+1)-XIOLD(I-1))
         SLOPE(I) = 0.0D0
C 
C  Evaluate function to be equidistributed
C 
         DO 300 J = 1, NTOL
            JJ = JTOL(J)
            JZ = LTOL(J)
            SLOPE(I) = MAX(SLOPE(I),(ABS(UHIGH(JJ,2)-UHIGH(JJ,1))
     *                 *WGTMSH(J)*ONEOVH/(1.0D0+ABS(Z(JZ,I))))**ROOT(J))
  300    CONTINUE
C 
C  Accumulate approximate integral of function to be equidistributed.
C 
         TEMP = SLOPE(I)*HIOLD
         SLPHMX = MAX(SLPHMX,TEMP)
         ACCUM(I+1) = ACCUM(I) + TEMP
         IFLIP = -IFLIP
  320 CONTINUE
C 
      AVRG = ACCUM(NOLD+1)/DBLE(NOLD)
      DEGEQU = AVRG/MAX(SLPHMX,50.0D0*X02AJF())
C 
C  NACCUM = expected N to achieve 0.1 x (user requested tolerances).
C 
      NACCUM = ACCUM(NOLD+1) + 1.0D0
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a,e12.5,a)')
     *     'DEB Degree of equidistribution is ', DEGEQU,
     *     ' for current mesh'
         CALL X04BAF(IOUT,REC)
         WRITE (REC,FMT='(a,i6)')
     *     'DEB Predicted number of meshpoints required ', N + 1
         CALL X04BAF(IOUT,REC)
      END IF
C 
C  Decide if mesh selection is worthwhile (otherwise, halve).
C 
      IF (AVRG.LT.50.0D0*X02AJF() .OR. DEGEQU.GE.0.5D0) GO TO 160
C 
C  NMX assures mesh has at least half as many subintervals as the
C  previous mesh.
C 
      NMX = MAX(NOLD+1,NACCUM)/2
C 
C  This assures that halving will be possible later (for error
C  estimate).
C 
      NMAX2 = NMAX/2
C 
C  The mesh is at most halved.
C 
      N = MIN(NMAX2,NOLD,NMX)
  340 CONTINUE
      NOLDP1 = NOLD + 1
      IF (N.LT.NFXP1) N = NFXP1
      MSHINF(NUMBER) = MSHINF(NUMBER) + 1
C 
C  If the new mesh is smaller than the old mesh set MSHINF(NUMBER)
C  so that the next call to  NEWMSH  will produce a halved
C  mesh. If N .eq. NOLD / 2 increment MSHINF(ALTER) so there can not
C  be an infinite loop alternating between N and N/2 points.
C 
      IF (N.LT.NOLD) MSHINF(NUMBER) = MSHINF(MLIMIT)
      IF (N.GT.NOLD/2) MSHINF(ALTER) = 1
      IF (N.EQ.NOLD/2) MSHINF(ALTER) = MSHINF(ALTER) + 1
      MSHINF(FLAG) = 0
C 
C  Having decided to generate a new mesh with N subintervals we now
C  do so, taking into account that the NFXPNT points in the array
C  FIXPNT must be included in the new mesh.
C 
      IN = 1
      ACCL = 0.0D0
      LOLD = 2
      XI(1) = ALEFT
      XI(N+1) = ARIGHT
      DO 460 I = 1, NFXP1
         IF (I.EQ.NFXP1) THEN
            ACCR = ACCUM(NOLDP1)
            LNEW = NOLDP1
            NREGN = N - IN
         ELSE
            DO 360 J = LOLD, NOLDP1
               LNEW = J
               IF (FIXPNT(I).LE.XIOLD(J)) GO TO 380
  360       CONTINUE
  380       CONTINUE
            ACCR = ACCUM(LNEW) + (FIXPNT(I)-XIOLD(LNEW))*SLOPE(LNEW-1)
            NREGN = (ACCR-ACCL)/ACCUM(NOLDP1)*DBLE(N) - 0.5D0
            NREGN = MIN(NREGN,N-IN-NFXP1+I)
            XI(IN+NREGN+1) = FIXPNT(I)
         END IF
         IF (NREGN.GT.0) THEN
            TEMP = ACCL
            TSUM = (ACCR-ACCL)/DBLE(NREGN+1)
            DO 440 J = 1, NREGN
               IN = IN + 1
               TEMP = TEMP + TSUM
               DO 400 L = LOLD, LNEW
                  LCARRY = L
                  IF (TEMP.LE.ACCUM(L)) GO TO 420
  400          CONTINUE
  420          CONTINUE
               LOLD = LCARRY
               XI(IN) = XIOLD(LOLD-1) + (TEMP-ACCUM(LOLD-1))
     *                  /SLOPE(LOLD-1)
  440       CONTINUE
         END IF
         IN = IN + 1
         ACCL = ACCR
         LOLD = LNEW
  460 CONTINUE
      MODE = 1
  480 CONTINUE
      NP1 = N + 1
      IF (IPRINT.LT.1) THEN
         WRITE (REC,FMT='(a,i6,a)') 'DEB New mesh of ', NP1, ' points:'
         CALL X04BAF(IOUT,REC)
         DO 500 I = 1, NP1, 5
            WRITE (REC,FMT='(a,5e14.6)') 'DEB ',
     *        (XI(J),J=I,MIN(I+4,NP1))
            CALL X04BAF(IOUT,REC)
  500    CONTINUE
      END IF
C 
      RETURN
      END
