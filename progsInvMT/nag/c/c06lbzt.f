      SUBROUTINE C06LBZ(CFUN,ZETA,RCIRC,EPREQ,EPMACH,NMAX,NCODE,EPEST,
     *                  NTCOF,TCOF,PARAMS,FENTRY)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     C06LBZ evaluates the normalized Taylor coefficients of a real
C     analytic function:
C
C     tcof(j+1) = (rcirc**j) * (j-th derivative of cfun(z) at z = zeta)
C                 divided by factorial(j)
C
C     for j = 0,1,2,3....nmax-1, to a uniform absolute accuracy epest
C     using function values of cfun(z) at points in the complex plane
C     lying on the circle of radius rcirc with center at z = zeta.
C
C     rcirc must be smaller than the radius of convergence of the Taylor
C     series. The problem has to be reformulated should cfun(z) happen
C     to be an odd function of (z - zeta), that is, if the relation
C
C        -cfun(-(z-zeta))=cfun(z-zeta)
C
C     is an identity.
C
C     C06LBZ is derived from the subroutine ENTCRE in the package WEEKS
C     by B.S. Garbow, G. Giunta, J.N. Lyness and A. Murli, Algorithm
C     662: A Fortran software package for the numerical inversion of the
C     Laplace Transform based on Weeks' method, ACM Trans. Math.
C     Software, 14, pp 171-176 (1988). This is a minor modification of
C     the subroutine ENTCRE by J.N. Lyness and G. Sande, Algorithm 413:
C     ENTCAF and ENTCRE, Comm. ACM, 14, pp 669-675 (1971). ENTCRE is a
C     special version of ENTCAF for use when zeta is real and also
C     cfun(z) is real when z is real.
C
C     C06LBZ call the NAG auxiliary routines C06EBV (for one-pass of a
C     radix-2 Hermitian-to-real FFT) and C06EAX (for a bit-reversal
C     permutation) instead of the routine HFTCOF. As a result the arrays
C     work and sintab are not needed.
C
C     **input parameters**
C
C     cfun     name of complex function subprogram.
C     zeta     real point about which taylor expansion is required.
C     rcirc    radius (real).
C     epreq    the absolute accuracy (real) to which the normalised
C              taylor coefficients, tcof(j), are required.
C     epmach   the machine accuracy parameter (real) (or an upper bound
C              on the relative accuracy of quantities likely to be
C              encountered).
C     nmax     physical upper limit on the size and length of the
C              calculation; the maximum number of coefficients
C              calculated will be that power of two less than or equal
C              to nmax. nmax is assumed to be at least 8.
C     ncode    .ge.0  the routine will do as well as  it can.
C              .lt.0  the routine will abort at an early stage if the
C                     required accuracy cannot be attained because of
C                     round-off error.
C     params   real parameter vector passed to cfun. (see note(2)
C              below.)
C     fentry   entry point passed to cfun. (see note(2) below.)
C
C     ** output parameters **
C
C     ncode    result status indicator. takes one of five values as
C              follows,
C                = +1. converged normally.
C                = -1. did not converge; no round-off error trouble.
C                = +2. converged, but with a higher tolerance set by
C                      the round-off level. ( epest.gt.epreq )
C                = -2. did not converge in spite of higher tolerance
C                      set by round-off level.
C                =  0. run was aborted because epreq is unattainable
C                      due to round-off level and input ncode is
C                      negative.
C     epest     estimate of actual uniform absolute accuracy in all
C               tcof (except, if ncode.eq.0, estimate of round-off
C               level).
C     ntcof     number of nontrivial values of  tcof  actually
C               calculated; they are based on  ntcof/2+2  calls of cfun
C               (three calls were for purely real argument).
C     tcof      approximations to the normalised taylor coefficients,
C               except when output ncode = 0.
C
C     note(1)   ncode is used both as input and output parameter;
C               normally it retains the value  +1  and need not be
C               reset between normal runs.
C     note(2)   cfun(z,params,fentry) is a user provided complex valued
C               function subprogram with a complex valued argument.
C               Possible parameters for cfun and/or callable
C               subprogram from cfun are communicated through params and
C               fentry.
C
C     ** bookkeeping parameters for stage one **
C
C     nconv   1  convergence achieved.
C            -1  no convergence achieved.
C     nround  1  no round off trouble observed.
C             2  round off trouble observed.
C     nabort  0  update tolerance and continue on appearance of round
C                off trouble.
C             1  terminate when round off trouble observed.
C     exact   the exact value of tcof(1) which is cfun(zeta).
C     safety  this is a safety factor by which the routine avoids the
C             round off level; it is set to 10.0 and appears only in
C             the combination (safety*epmach). To alter this factor, or
C             to remove the round off error guard completely, the user
C             need only adjust the input parameter epmach appropriately.
C
C     ** quantities calculated in stage three(a) **
C
C     this is the first part of iteration number ntcof.
C     this stage is now empty.
C
C     ** quantities calculated in stage three(b) **
C
C     iterations are numbered 8,16,32... at the end of iteration number
C     ntcof, the ntcof/2 + 1 complex function values at abscissae
C     regularly spaced on upper half of circle are stored in the tcof
C     vector as follows:
C
C      tcof(j+1)       =    real   part of cfun(z(j))
C                                                 j=0,1,2,....ntcof/2.
C      tcof(ntcof-j+1) = imaginary part of cfun(z(j))
C                                                 j=1,2,...(ntcof/2-1).
C     where
C      z(j)  =  zeta + rcirc*cexp(2*pi*eye*j/ntcof)
C     this involves  a rearrangement of the ntcof/4 + 1 function values
C     available at the start of the iteration and the calculation of a
C     further ntcof/4 function values; in addition fmax and approx are
C     calculated. These are
C      fmax     maximum modulus of the function values so far
C               encountered.
C      approx   an approximation to tcof(1) based on these function
C               values.
C
C     ** quantities calculated at stage three(c) **
C
C     error1  current value of the error = abs(approx-exact).
C     error2,error3,error4  values of error at end of three previous
C             iterations.
C     epmach  machine accuracy parameter. (input parameter)
C     epreq   required accuracy. (input parameter)
C     epro    highest accuracy reasonably attainable in view of the size
C             of the function values so far encountered.
C             (=10.0*epmach*fmax)
C     epcof   currently required accuracy (=amax1(epreq,epro)).
C     epest   estimate of current accuracy. (the maximum of epro and a
C             function of errors 1,2,3 and 4 ) (output parameter)
C
C     ** convergence and termination checks in stage three(c) **
C
C     (1)  uses fmax to raise epcof above round off level; if this is
C          necessary and the input value of ncode is negative, it
C          terminates setting ncode = 0.
C     (2)  uses approx to evaluate convergence of tcof(1) towards exact;
C          it may assign convergence and go to stage four(a) setting
C          ncode = +1 or +2.
C     (3)  uses nmax to check physical limit; if this has been reached,
C          it goes to stage four(a) setting ncode = -1 or -2.
C     (4)  otherwise continues next iteration by going to stage three.
C
C     **  calculation of first ntcof taylor coefficients in stage
C                                                             four(a) **
C
C     The fft calculates the neccessary summations except for dividing
C     by ntcof.
C
C     **  setting of remaining taylor coefficients in stage four(b)  **
C
C     the convergence criterion allows us to infer that the normalized
C     taylor coefficients of order greater than  ntcof  are zero to
C     accuracy epest. They are evaluated as being exactly zero.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPEST, EPMACH, EPREQ, RCIRC, ZETA
      INTEGER           NCODE, NMAX, NTCOF
C     .. Array Arguments ..
      DOUBLE PRECISION  PARAMS(3), TCOF(NMAX)
C     .. Function Arguments ..
      COMPLEX*16        CFUN, FENTRY
      EXTERNAL          CFUN, FENTRY
C     .. Local Scalars ..
      COMPLEX*16        FVAL, ZVAL
      DOUBLE PRECISION  APPROX, EP32, EP42, EPCOF, EPMIN, EPRO, ERROR1,
     *                  ERROR2, ERROR3, ERROR4, EXACT, FMAX, FVALIM,
     *                  FVALRE, RCOS, RSIN, SAFETY, SCALE, SUPPER, TPN,
     *                  TWOPI
      INTEGER           ITWO, J, JCONJ, JRCONJ, JREFL, NABORT, NCONV,
     *                  NDISP, NPREV, NROUND
C     .. Local Arrays ..
      INTEGER           TWOS(21)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. External Subroutines ..
      EXTERNAL          C06EAX, C06EBV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, MAX, MIN, DBLE, SIN
C     .. Executable Statements ..
C
C     ***   stage one   ***
C
C     initialise bookkeeping parameters and exact value of tcof(1).
C
      TWOPI = 2.0D0*X01AAF(ZERO)
      NROUND = 1
      NABORT = 0
      IF (NCODE.LT.0) NABORT = 1
      EPCOF = EPREQ
      SAFETY = 10.0D0
      ZVAL = DCMPLX(ZETA,ZERO)
      FVAL = CFUN(ZVAL,PARAMS,FENTRY)
      FVALRE = DBLE(FVAL)
      EXACT = FVALRE
C
C     ***   stage two   ***
C
C     first three iterations ( those with ntcof = 1,2,4 ).
C
      ZVAL = DCMPLX(ZETA+RCIRC,ZERO)
      FVAL = CFUN(ZVAL,PARAMS,FENTRY)
      FVALRE = DBLE(FVAL)
      APPROX = FVALRE
      FMAX = ABS(FVALRE)
      TCOF(1) = FVALRE
      ERROR3 = ABS(APPROX-EXACT)
      ZVAL = DCMPLX(ZETA-RCIRC,ZERO)
      FVAL = CFUN(ZVAL,PARAMS,FENTRY)
      FVALRE = DBLE(FVAL)
      APPROX = 0.5D0*(APPROX+FVALRE)
      FMAX = MAX(FMAX,ABS(FVALRE))
      TCOF(3) = FVALRE
      ERROR2 = ABS(APPROX-EXACT)
      ZVAL = DCMPLX(ZETA,RCIRC)
      FVAL = CFUN(ZVAL,PARAMS,FENTRY)
      FVALRE = DBLE(FVAL)
      FVALIM = DIMAG(FVAL)
      APPROX = 0.5D0*(APPROX+FVALRE)
      FMAX = MAX(FMAX,ABS(FVAL))
      TCOF(2) = FVALRE
      TCOF(4) = FVALIM
      ERROR1 = ABS(APPROX-EXACT)
      NTCOF = 4
      EPRO = FMAX*SAFETY*EPMACH
      IF (EPRO.LT.EPCOF) GO TO 20
      EPCOF = EPRO
      NROUND = 2
      IF (NABORT.EQ.0) GO TO 20
      NCODE = 0
      EPEST = EPRO
      GO TO 200
C
C     ***   stage three   ***
C
C     commence iteration number ntcof.
C
   20 CONTINUE
      NPREV = NTCOF
      NTCOF = 2*NTCOF
C
C     ***   stage three(a)   ***
C
C     this stage is now empty
C
C     ***   stage three(b)   ***
C
C     update list of function values in tcof, calculate fmax and approx.
C
      DO 40 J = NPREV - 1, 1, -1
         TCOF(2*J+1) = TCOF(J+1)
   40 CONTINUE
      SUPPER = ZERO
      TPN = TWOPI/DBLE(NTCOF)
      DO 60 J = 1, (NPREV/2) - 1, 2
         RSIN = RCIRC*SIN(TPN*J)
         RCOS = RCIRC*COS(TPN*J)
         JCONJ = NTCOF - J
         ZVAL = DCMPLX(ZETA+RCOS,RSIN)
         FVAL = CFUN(ZVAL,PARAMS,FENTRY)
         FVALRE = DBLE(FVAL)
         FVALIM = DIMAG(FVAL)
         SUPPER = SUPPER + FVALRE
         FMAX = MAX(FMAX,ABS(FVAL))
         TCOF(J+1) = FVALRE
         TCOF(JCONJ+1) = FVALIM
         JREFL = NPREV - J
         JRCONJ = NTCOF - JREFL
         ZVAL = DCMPLX(ZETA-RCOS,RSIN)
         FVAL = CFUN(ZVAL,PARAMS,FENTRY)
         FVALRE = DBLE(FVAL)
         FVALIM = DIMAG(FVAL)
         SUPPER = SUPPER + FVALRE
         FMAX = MAX(FMAX,ABS(FVAL))
         TCOF(JREFL+1) = FVALRE
         TCOF(JRCONJ+1) = FVALIM
   60 CONTINUE
      APPROX = 0.5D0*APPROX + SUPPER/NPREV
C
C     ***   stage three(c)   ***
C
C     convergence and termination check.
C
      ERROR4 = ERROR3
      ERROR3 = ERROR2
      ERROR2 = ERROR1
      ERROR1 = ABS(APPROX-EXACT)
      EPRO = FMAX*SAFETY*EPMACH
      IF (EPRO.LT.EPCOF) GO TO 80
      EPCOF = EPRO
      NROUND = 2
      IF (NABORT.EQ.0) GO TO 80
      NCODE = 0
      EPEST = EPRO
      GO TO 200
   80 CONTINUE
      ERROR4 = MAX(ERROR4,EPRO)
      ERROR3 = MAX(ERROR3,EPRO)
      EP42 = ERROR2*((ERROR2/ERROR4)**(4.0D0/3.0D0))
      EP32 = ERROR2*((ERROR2/ERROR3)**2)
      EPMIN = MIN(ERROR2,EP32,EP42)
      EPEST = MAX(ERROR1,EPMIN,EPRO)
      IF (EPEST.GT.EPCOF) GO TO 100
      NCONV = 1
      GO TO 120
  100 CONTINUE
      IF (2*NTCOF.LE.NMAX) GO TO 20
      NCONV = -1
C
C     ***   stage four(a)   ***
C
C     calculation of first ntcof taylor coefficients using f.f.t.
C
  120 CONTINUE
      NCODE = NCONV*NROUND
C
C     call NAG routines to perform f.f.t.
C
      NDISP = NTCOF
      ITWO = 0
  140 CONTINUE
      NDISP = NDISP/2
      ITWO = ITWO + 1
      TWOS(ITWO) = 2
C
C     perform one pass of radix-2 Hermitian-to-real f.f.t.
C
      CALL C06EBV(TCOF,NTCOF,TCOF(NDISP+1),NTCOF-NDISP,NDISP)
      IF (NDISP.GT.1) GO TO 140
C
C     re-order results
C
      TWOS(ITWO+1) = 0
      CALL C06EAX(TCOF,NTCOF,TWOS)
      SCALE = ONE/NTCOF
      DO 160 J = 1, NTCOF
         TCOF(J) = TCOF(J)*SCALE
  160 CONTINUE
C
C     ***   stage four(b)   ***
C
C     setting of remaining taylor coefficients.
C
      DO 180 J = NTCOF + 1, NMAX
         TCOF(J) = ZERO
  180 CONTINUE
  200 CONTINUE
      RETURN
C     end of C06LBZ
      END
