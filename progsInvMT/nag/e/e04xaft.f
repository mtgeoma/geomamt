      SUBROUTINE E04XAF(MSGLVL,N,EPSRF,X,MODE,OBJFUN,LHES,HFORW,FX,GRAD,
     *                  HCNTRL,HESIAN,IWARN,WORK,IUSER,USER,INFO,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13B REVISED. IER-659 (AUG 1988).
C     MARK 14 REVISED. IER-725 (DEC 1989).
C     MARK 16 REVISED. IER-1106 (JUL 1993).
C***********************************************************************
C     E04XAF (FDCALC)  computes a set of forward-difference intervals
C     for the n-variable function objfun at the point x by repeatedly
C     calling subroutine  E04XAZ (CHCORE).
C
C     Full documentation of  fdcalc  is contained in report SOL 83-6,
C     Documentation of FDCORE and FDCALC, by P.E. Gill, W. Murray,
C     M.A. Saunders and M.H. Wright, Department of Operations Research,
C     Stanford University, Stanford, California 94305, June 1983.
C
C     The input parameters of E04XAF are ...
C
C     msglvl   integer
C              -Controls the level of printout.
C               msglvl = 0  No printout.
C               msglvl = 1  Summary printout for each variable and
C                           warning messages.
C               msglvl = 2  Full debug printout from E04XAF and CHCORE.
C
C     n        integer
C              -The number of variables.
C
C     epsrf    real
C              -A good bound on the relative error in computing fx at x.
C               If epsrf is negative on entry or if the users value of
C               epsrf is too small or too large then the default value
C               of e**0.9 is used, where e is the machine precision.
C
C     x        real array of dimension at least (n).
C              -The point at which the intervals are to be computed.
C
C     mode     integer
C              -If mode = 0 E04XAF returns the gradient vector and
C                           hessian diagonal given function only.
C               If mode = 1 E04XAF returns the hessian matrix given
C                           function and gradients.
C               If mode = 2 E04XAF returns the gradient vector and
C                           hessian matrix given function only.
C
C     objfun   subroutine provided by the user.
C              -Objfun must calculate the objective function ( and
C               gradients if mode = 1 ). Objfun must be declared as
C               external in the routine that calls E04XAF.
C               Its specification is:
C               subroutine objfun( mode, n, x, fx grad, nstate)
C
C     lhes     integer
C              -The first dimension of array hesian as declared in the
C               calling program.
C               lhes .ge. n
C
C     the input/output parameters of E04XAF are ...
C
C     hforw    real array of dimension at least (n).
C              -Before entry, hforw should contain the initial trial
C               intervals for:-
C                  variable x(i)    if mode = 0 or mode = 2
C                  gradient grad(i) if mode = 1
C               If hforw(i) .le. 0.0 then the initial trial interval is
C               computed by E04XAF.
C               On exit hforw contains the forward-difference intervals.
C
C     the output parameters of E04XAF are ...
C
C     fx       real
C              -The value of the function at x.
C
C     grad     real array of dimension at least (n).
C              -The approximate gradient vector at x (computed by
C               forward-differences with  hforw if mode = 0 or 2).
C
C     hcntrl   real array of dimension at least (n).
C              -The central-difference intervals.
C
C     hesian   real array of dimension (lhes,n) if mode = 1 or mode = 2.
C                ''  ''   ''  ''       (lhes,1) if mode = 0.
C              -If mode = 1 or mode = 2 then the leading n by n part of
C               this array will contain the full hessian matrix.
C               If mode = 0 then the first n elements of this array will
C               contain the hessian diagonals.
C
C     iwarn    integer
C              -On successful exit iwarn = 0.
C               If the value of epsrf on entry is too small or too large
C               then iwarn is set to 1 or 2 respectively on exit.
C
C     the workspace parameters of E04XAF are ...
C
C     work     real array of dimension at least (n**2 + n) if mode = 1
C              mode = 2.
C              real array of dimension at least (n) if mode = 0.
C              -Used as workspace.
C
C     iuser    integer array of length at least (1)
C              -This array is not used by E04XAF, but is passed directly
C               to user supplied routine objfun and may be used to
C               supply information to objfun.
C
C     user     real array of length at least (1)
C              -This array is not used by E04XAF, but is passed directly
C               to user supplied routine objfun and may be used to
C               supply information to objfun.
C
C     the diagnostic parameters of E04XAF are ...
C
C     info     integer array of dimension at least (n)
C              -The i-th component indicates the status of the interval
C               for variable i (if mode = 0 or 2) or gradient i (if
C               mode = 1). info(i) = 0  means that the interval was
C               satisfactory.
C
C     ifail    integer
C              -Before entry ifail must be set to 0, -1 or 1 (see
C               Chapter p01). The recommended value is -1.
C               On successful exit ifail = 0.
C
C     Original Fortran 66 Version 2.1 of June 1983.  (SOL)
C     Fortran 77 version dated 17-May-1985.          (SOL)
C     Edited to output full hessian matrix. Mark 12. (NAG)
C
C***********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, THREE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0,THREE=3.0D+0)
      DOUBLE PRECISION  POINT9
      PARAMETER         (POINT9=0.9D+0)
      DOUBLE PRECISION  ETA, RHO
      PARAMETER         (ETA=1.0D+0,RHO=1.0D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSRF, FX
      INTEGER           IFAIL, IWARN, LHES, MODE, MSGLVL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  GRAD(N), HCNTRL(N), HESIAN(LHES,*), HFORW(N),
     *                  USER(*), WORK(*), X(N)
      INTEGER           INFO(N), IUSER(*)
C     .. Subroutine Arguments ..
      EXTERNAL          OBJFUN
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANUM, CDEST, DENOM, EPSA, EPSMCH, EPSPT9, EPSR,
     *                  ERRBND, F, F1, F2, FDEST, FEPS, FORG, FXHI,
     *                  FXHIHJ, FXHJ, H, HOPT, HPHI, SDEST, TWOTHI,
     *                  XISAVE, XJSAVE
      INTEGER           I, ICOUNT, INFORM, ITER, ITMAX, J, NADV, NFI,
     *                  NFSVE, NSTATE, NUMF
      LOGICAL           DEBUG, DONE, FIRST, HEADNG, IERR, OVRFLW
      CHARACTER*6       SRNAME
C     .. Local Arrays ..
      CHARACTER*120     REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  F06BLF
      INTEGER           P01ABF
      EXTERNAL          F06BLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04XAZ, X02ZAZ, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              SRNAME/'E04XAF'/
C     .. Executable Statements ..
C
C     Compute the machine-dependent constants.
C
      CALL X02ZAZ
C
C     Set NADV, the output unit for advisory messages
C
      NADV = WMACH(11)
      NOUT = WMACH(11)
      EPSMCH = WMACH(3)
      EPSPT9 = EPSMCH**POINT9
      TWOTHI = TWO/THREE
      EPSR = EPSRF
      HEADNG = .TRUE.
C
C     Initialise IERR to false, it will be set to true if an input
C     parameter is invalid or if any element on exit from E04XAZ has
C     a non-zero INFORM value.
C     First check input parameter values.
C
      IERR = .FALSE.
      IF (N.LT.1) THEN
         IERR = .TRUE.
         WRITE (REC,FMT=99999)
      ELSE IF (MODE.NE.0 .AND. MODE.NE.1 .AND. MODE.NE.2) THEN
         IERR = .TRUE.
         WRITE (REC,FMT=99994)
      ELSE IF (LHES.LT.N) THEN
         IERR = .TRUE.
         WRITE (REC,FMT=99993)
      END IF
      IF (IERR) THEN
         IFAIL = P01ABF(IFAIL,1,SRNAME,1,REC(1))
         RETURN
      END IF
C
      ITMAX = 3
C
C     Set the local logical variable  debug  to  .true.  if  msglvl = 2.
C
      IF (MSGLVL.NE.2) DEBUG = .FALSE.
      IF (MSGLVL.EQ.2) DEBUG = .TRUE.
C
C     Compute the function at the base point.
C     NUMF represents the number of function evaluations. NSTATE is set
C     to 1 on the first call of objfun and is 0 for all subsequent calls
C     this allows the user to perform certain calculations once only.
C
      NUMF = 0
      NSTATE = 1
      CALL OBJFUN(MODE,N,X,FX,GRAD,NSTATE,IUSER,USER)
      NSTATE = 0
      NUMF = NUMF + 1
      IF (MODE.LT.0) GO TO 180
C
C     Test the value of EPSR. If EPSR is too small or too large
C     output message.
C
      IWARN = 0
      IF (EPSR.LE.ZERO) THEN
         EPSR = EPSPT9
      ELSE
         IF (FX.GE.EPSMCH**TWOTHI) THEN
            FEPS = FX*EPSR
            IF (FEPS.GE.ONE) THEN
               IWARN = 2
               IF (MSGLVL.GE.1) THEN
                  WRITE (REC,FMT=99996)
                  CALL X04BAF(NADV,REC(1))
               END IF
               EPSR = EPSPT9
            END IF
            IF (FEPS.LT.EPSMCH) THEN
               IWARN = 1
               IF (MSGLVL.GE.1) THEN
                  WRITE (REC,FMT=99995)
                  CALL X04BAF(NADV,REC(1))
               END IF
               EPSR = EPSPT9
            END IF
         END IF
      END IF
C
C     If  debug  is  .true.,  a detailed intermediate printout will
C     be produced.
C
      EPSA = EPSR*(ONE+ABS(FX))
      IF (DEBUG) THEN
         WRITE (REC,FMT=99992) EPSR, FX
         CALL X04BAF(NADV,REC(1))
      END IF
C
C     ==================================================================
C     Loop over each of the components of  x.
C     ==================================================================
      DO 40 I = 1, N
C
C        Save this component of  X  and initialize  FIRST.
C        The parameter  FIRST  indicates which portion of E04XAZ is
C        to be executed.  It must be set to true before the first call
C        of  E04XAZ for each variable.
C
         XISAVE = X(I)
         FIRST = .TRUE.
         ITER = 0
         CDEST = ZERO
         SDEST = ZERO
         NFSVE = NUMF
C
         IF (DEBUG) THEN
            WRITE (REC,FMT=99991) I, XISAVE
            CALL X04BAF(NADV,REC(1))
         END IF
C
         IF (MODE.EQ.1) THEN
            FORG = GRAD(I)
            EPSA = EPSR*(ONE+ABS(GRAD(I)))
         ELSE
            FORG = FX
         END IF
C
C        Set the interval value HOPT, this is the interval that would
C        be optimal for a well scaled problem.
C
         IF (HFORW(I).GT.ZERO) THEN
            HOPT = HFORW(I)
            H = HOPT
            IF (MODE.EQ.0) THEN
               H = RHO*HOPT
            END IF
         ELSE
            IF (MODE.EQ.0) THEN
               HOPT = TWO*(ETA+ABS(XISAVE))*SQRT(EPSR)
               H = RHO*HOPT
            ELSE IF (MODE.EQ.1) THEN
               HOPT = TWO*(ETA+ABS(XISAVE))*SQRT(EPSR)
               H = HOPT
            ELSE IF (MODE.EQ.2) THEN
               HOPT = TWO*(ETA+ABS(XISAVE))*((EPSR)**0.25D+0)
               H = HOPT
            END IF
         END IF
C
C        Compute the next trial interval.
C        repeat
C
   20    CONTINUE
C
C        Compute the function at the new trial points.
C
         X(I) = XISAVE + H
         CALL OBJFUN(MODE,N,X,F1,WORK,NSTATE,IUSER,USER)
         NUMF = NUMF + 1
         IF (MODE.LT.0) THEN
            GO TO 180
         ELSE IF (MODE.EQ.1) THEN
            F1 = WORK(I)
         END IF
         X(I) = XISAVE + H + H
         CALL OBJFUN(MODE,N,X,F2,WORK,NSTATE,IUSER,USER)
         NUMF = NUMF + 1
         IF (MODE.LT.0) THEN
            GO TO 180
         ELSE IF (MODE.EQ.1) THEN
            F2 = WORK(I)
         END IF
C
         CALL E04XAZ(DEBUG,DONE,FIRST,EPSA,EPSR,FORG,INFORM,ITER,ITMAX,
     *               CDEST,FDEST,SDEST,ERRBND,F1,F2,H,HOPT,HPHI)
C
C
C        until  DONE
         IF ( .NOT. DONE) GO TO 20
C
C        ===========================================================
C        Exit for this variable.
C        Assign the approximate gradient and compute the Hessian as
C        required.
C        ===========================================================
C
         INFO(I) = INFORM
         IF (INFO(I).GT.0) IERR = .TRUE.
C
         IF (MODE.EQ.0 .OR. MODE.EQ.2) THEN
            GRAD(I) = CDEST
            HESIAN(I,1) = SDEST
         ELSE
            HESIAN(I,1) = CDEST
         END IF
         HFORW(I) = HOPT
         HCNTRL(I) = HPHI
         X(I) = XISAVE
         NFI = NUMF - NFSVE
C
C        A positive value of  msglvl  indicates that a summary printout
C        for each variable will be given by  E04XAF.
C
         IF (MSGLVL.GT.0) THEN
            IF (HEADNG) THEN
               IF (MODE.EQ.1) THEN
                  WRITE (REC,FMT=99988)
                  CALL X04BAY(NADV,5,REC)
               ELSE
                  WRITE (REC,FMT=99990)
                  CALL X04BAY(NADV,5,REC)
               END IF
               HEADNG = .FALSE.
            END IF
            IF (MODE.EQ.1) THEN
               WRITE (REC,FMT=99987) I, X(I), HFORW(I), HCNTRL(I),
     *           ERRBND, HESIAN(I,1), NFI, INFO(I)
               CALL X04BAF(NADV,REC(1))
            ELSE
               WRITE (REC,FMT=99989) I, X(I), HFORW(I), HCNTRL(I),
     *           ERRBND, GRAD(I), HESIAN(I,1), NFI, INFO(I)
               CALL X04BAF(NADV,REC(1))
            END IF
         END IF
C
   40 CONTINUE
C
C     Compute the approx to the hessian given function and gradients
C     Reference: Practical Optimization by Gill,P.E., Murray,W. and
C     Wright,M.H. pages 54 - 56.
C     Forward-difference approx to the j-th column of hessian is given
C     by y(j) = ( 1.0/h(j) )*( g(x + h(j)*e(j)) - g(x) )
C
      IF (MODE.EQ.1) THEN
         ICOUNT = 0
         DO 80 I = 1, N
            XISAVE = X(I)
            X(I) = X(I) + HFORW(I)
            NUMF = NUMF + 1
            CALL OBJFUN(MODE,N,X,F,WORK,NSTATE,IUSER,USER)
            DO 60 J = 1, N
               ICOUNT = ICOUNT + 1
               WORK(J) = (WORK(J)-GRAD(J))/HFORW(I)
               WORK(N+ICOUNT) = WORK(J)
   60       CONTINUE
            X(I) = XISAVE
   80    CONTINUE
C
C        To make hessian symmetric perform  G = 0.5*( Y + trans(Y) )
C
         DO 120 J = 1, N
            DO 100 I = J, N
               HESIAN(I,J) = 0.5D0*(WORK(N+(J-1)*N+I)+WORK(N+(I-1)*N+J))
               IF (I.NE.J) THEN
                  HESIAN(J,I) = HESIAN(I,J)
               END IF
  100       CONTINUE
  120    CONTINUE
C
      ELSE IF (MODE.EQ.2) THEN
C
C        Compute hessian given function only.
C        Reference: Practical Optimization by Gill,P.E., Murray,W. and
C        Wright,M.H. pages 54 - 56.
C
C
         DO 160 I = 1, N
            DO 140 J = 1, N
               IF (I.LE.J) THEN
                  XISAVE = X(I)
                  XJSAVE = X(J)
                  X(I) = X(I) + HCNTRL(I)
                  X(J) = X(J) + HCNTRL(J)
                  NUMF = NUMF + 1
                  CALL OBJFUN(MODE,N,X,FXHIHJ,WORK,NSTATE,IUSER,USER)
                  X(I) = XISAVE
                  X(J) = XJSAVE + HCNTRL(J)
                  NUMF = NUMF + 1
                  CALL OBJFUN(MODE,N,X,FXHJ,WORK,NSTATE,IUSER,USER)
                  X(J) = XJSAVE
                  X(I) = XISAVE + HCNTRL(I)
                  NUMF = NUMF + 1
                  CALL OBJFUN(MODE,N,X,FXHI,WORK,NSTATE,IUSER,USER)
                  X(I) = XISAVE
C
                  ANUM = FXHIHJ - FXHJ - FXHI + FX
                  DENOM = HCNTRL(I)*HCNTRL(J)
                  HESIAN(I,J) = F06BLF(ANUM,DENOM,OVRFLW)
               ELSE
                  HESIAN(I,J) = HESIAN(J,I)
               END IF
  140       CONTINUE
  160    CONTINUE
      END IF
C
C     If one or more variables have a nonzero inform value
C     set ifail to 2
C
      IF (IERR) THEN
         WRITE (REC,FMT=99998)
         IFAIL = P01ABF(IFAIL,2,SRNAME,1,REC(1))
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
C     The user has requested termination by setting mode negative in
C     objfun.
C
  180 CONTINUE
      WRITE (REC,FMT=99997)
      IFAIL = P01ABF(IFAIL,MODE,SRNAME,1,REC(1))
      RETURN
C
C     End of E04XAF (FDCALC).
C
99999 FORMAT (' Error return from E04XAF because N is less than 1')
99998 FORMAT (' One or more variables have a nonzero INFO value')
99997 FORMAT (' User requested termination by setting MODE negative in',
     *       ' routine OBJFUN.')
99996 FORMAT (' WARNING EPSRF is too large and has been set to default',
     *       ' value.')
99995 FORMAT (' WARNING EPSRF is too small and has been set to default',
     *       ' value')
99994 FORMAT (' Error return from E04XAF because MODE is out of range')
99993 FORMAT (' Error return from E04XAF because LHES is less than N')
99992 FORMAT (' //E04XAF//  Epsrf and initial function',1P,2D16.6)
99991 FORMAT (' //E04XAF//  Variable',I5,' = ',1P,D16.6)
99990 FORMAT (//' Output from routine E04XAF.',/'    J      X(J)   F. ',
     *       'dif. int.   C. dif. int.     Error est.     Grad. est. H',
     *       'ess diag est. fun evals.  info(j)',/)
99989 FORMAT (I5,1P,D10.2,1P,5D15.6,I9,I8)
99988 FORMAT (//' Output from routine E04XAF.',/'    J      X(J)   F. ',
     *       'dif. int.   C. dif. int.     Error est. Hess diag est. f',
     *       'un evals.  info(j)',/)
99987 FORMAT (I5,1P,D10.2,1P,4D15.6,I9,I8)
      END
