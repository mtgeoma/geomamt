      SUBROUTINE D02PDY(TNOW,Y,YP,TSTG,YSTG,YPSTG,HTRY,WEIGHT,CUTBAK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE STEPA $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      To calculate an "on-scale" step size for phase 2 of
C                the initial step size computation.
C
C  Input:        TNOW, Y(*), YP(*), TSTG, YSTG(*), YPSTG(*)
C  Input/output: HTRY, WEIGHT
C  Output:       CUTBAK
C
C  Common:       Initializes:    none
C                Reads:          /AD02PD/ TND, NEQ
C                                /ED02PD/ STBRAD, RS1, RS4
C                                /GD02PD/ RNDOFF
C                Alters:         none
C
C  Comments:
C  =========
C  This subroutine is used during the first three stages of the first
C  step. A Lipschitz constant L for the differential equation in
C  autonomous form is approximated, and the product abs(HTRY)*L is
C  compared to an approximate radius, STBRAD, of the stability region
C  of the method. The step size is reduced as necessary, within a range
C  specified by the step size control parameters RS1 and RS4, to assure
C  stability and give some confidence in the error estimator. If HTRY is
C  reduced, CUTBAK is set .TRUE..
C
C  Y(*) and YP(*) contain the solution and its derivative at TNOW and
C  similarly YSTG(*) and YPSTG(*) contain approximations at TSTG.
C
C  Normally the weights used in the control of the error depend on the
C  size of the solution at the beginning and at the end of the step, but
C  at this time we do not have a solution at the end of the step.  Each
C  stage YSTG(*) of the Runge - Kutta process represents a low order
C  approximation to the solution at TSTG.  Because the initial value of
C  WEIGHT(*) provided in the first phase of the scheme is based only on
C  the solution at T and THRES(*), it is continually updated in D02PDY
C  to account for the size of the solution throughout the step as
C  revealed by the intermediate stages YSTG(*). Inside this subroutine
C  only, the differential equation is converted to autonomous form.
C  After the conversion, the end of the interval of integration, TND,
C  is used to define a suitable weight for the independent variable.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HTRY, TNOW, TSTG
      LOGICAL           CUTBAK
C     .. Array Arguments ..
      DOUBLE PRECISION  WEIGHT(*), Y(*), YP(*), YPSTG(*), YSTG(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, CUBRMC, DIR, DWARF, EXPON, HSTRT, MCHEPS,
     *                  RNDOFF, RS, RS1, RS2, RS3, RS4, SAFETY, SQRRMC,
     *                  STBRAD, TANANG, TINY, TND, TOLR, TOOSML, TSTRT
      INTEGER           LSTSTG, MAXTRY, NEQN, NSEC, ORDER, OUTCH
      LOGICAL           FSAL
C     .. Local Scalars ..
      DOUBLE PRECISION  ARGDIF, FDIFF, SCL, TDIFF, TWT, WT, YNRM, YSTGNM
      INTEGER           L
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /ED02PD/TOOSML, COST, SAFETY, EXPON, STBRAD,
     *                  TANANG, RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG,
     *                  MAXTRY, NSEC, FSAL
      COMMON            /GD02PD/MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC,
     *                  TINY, OUTCH
C     .. Save statement ..
      SAVE              /AD02PD/, /ED02PD/, /GD02PD/
C     .. Executable Statements ..
C
C  Update the weights to account for the current intermediate solution
C  approximation YSTG(*). Compute the sizes of Y(*) and YSTG(*) in the
C  new norm. The size of the Lipschitz constant is assessed by a
C  difference in the arguments Y(*), YSTG(*) and a difference in the
C  function evaluated at these arguments.
C
      YNRM = ZERO
      YSTGNM = ZERO
      ARGDIF = ZERO
      FDIFF = ZERO
      DO 20 L = 1, NEQN
         WT = MAX(WEIGHT(L),ABS(YSTG(L)))
         WEIGHT(L) = WT
         YNRM = MAX(YNRM,ABS(Y(L))/WT)
         YSTGNM = MAX(YSTGNM,ABS(YSTG(L))/WT)
         ARGDIF = MAX(ARGDIF,ABS(YSTG(L)-Y(L))/WT)
         FDIFF = MAX(FDIFF,ABS(YPSTG(L)-YP(L))/WT)
   20 CONTINUE
C
C  The transformation of the equation to autonomous form is done
C  implicitly.  The difference of the arguments must take into account
C  the difference between the values of the independent variable T and
C  TSTG. The difference of the corresponding component of the function
C  is zero because of the way the standard transformation is done.
C
      TDIFF = TSTG - TNOW
      TWT = ABS(TND-TNOW)
      YNRM = MAX(YNRM,ABS(TNOW)/TWT)
      YSTGNM = MAX(YSTGNM,ABS(TSTG)/TWT)
      ARGDIF = MAX(ARGDIF,ABS(TDIFF)/TWT)
C
C  The ratio FDIFF/ARGDIF is a lower bound for, and an approximation to,
C  a Lipschitz constant L for the differential equation written in
C  autonomous form. First we must ask if the difference ARGDIF is
C  significant in the precision available. If it appears to be, we
C  insist that abs(HTRY)*L be less than an approximate radius, STBRAD,
C  of the stability region of the method. This is more stringent than
C  necessary for stability, possibly a lot more stringent, but the aim
C  is to get an HTRY small enough that the error estimate for the step
C  is credible.  The reduction is required to be at least as much as the
C  step control parameter RS1. It is necessary to limit the reduction of
C  HTRY at any one time because we may be misled in the size of the
C  reduction that is appropriate due to nonlinearity of the differential
C  equation and to inaccurate weights caused by HTRY much too large. The
C  reduction is not permitted to be more than the step control parameter
C  RS4.
C
      CUTBAK = .FALSE.
      IF (ARGDIF.GT.RNDOFF*MAX(YNRM,YSTGNM)) THEN
         IF ((ABS(HTRY)*FDIFF).GT.(STBRAD*ARGDIF)) THEN
            SCL = (STBRAD*ARGDIF)/(ABS(HTRY)*FDIFF)
            SCL = MIN(SCL,RS1)
            SCL = MAX(SCL,RS4)
            HTRY = SCL*HTRY
            CUTBAK = .TRUE.
         END IF
      END IF
C
      RETURN
      END
