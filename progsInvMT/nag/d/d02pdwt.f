      SUBROUTINE D02PDW(F,NEQ,Y,TOL,WEIGHT,ZY,ZYP,ZERROR,ZYNEW,ZERRES,
     *                  ZSTAGE,IER)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE TRUERR $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:     Compute a running RMS measure of the true (global) error
C               for a general Runge-Kutta pair.
C
C
C  Input:        NEQ, Y(*), TOL, WEIGHT(*),
C  Input/output: ZY(*), ZYP(*), ZERROR(*)
C  Workspace:    ZYNEW(*), ZERRES(*), ZSTAGE(NEQ,*)
C  Output:       IER
C  External:     F
C
C  Common:       Initializes:    none
C                Reads:          /BD02PD/ T, HOLD
C                                /ED02PD/ TOOSML, ORDER, NSEC
C                                /GD02PD/ TINY
C                Alters:         /FD02PD/ MAXERR, LOCMAX, GNFCN
C
C  Comments:
C  =========
C  A secondary integration is performed using a fraction of the step
C  size of the primary integration. ZY(*) and ZYP(*) are the
C  approximate solution and first derivative of this secondary
C  integration. ZERRES(*) contains the error estimates for the
C  secondary integration. ZYNEW(*) and ZSTAGE(*,*) are workspace for
C  taking a step. The error assessment is computed using the difference
C  of the primary and secondary solutions at the primary integration
C  points as an estimate of the true error there.  The weights used are
C  those of the error test of the primary integration. This error
C  assessment is maintained in the vector ZERROR(*). MAXERR and LOCMAX
C  contain the maximum contribution to the assessment and its location,
C  respectively. The number of calls to F is counted by GNFCN.
C
C     .. Parameters ..
      DOUBLE PRECISION  PT1, TEN, DUMMY
      PARAMETER         (PT1=0.1D0,TEN=10.0D0,DUMMY=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IER, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  WEIGHT(*), Y(*), ZERRES(*), ZERROR(*),
     *                  ZSTAGE(NEQ,*), ZY(*), ZYNEW(*), ZYP(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, CUBRMC, DWARF, EXPON, H, HOLD, LOCMAX,
     *                  MAXERR, MCHEPS, RNDOFF, RS, RS1, RS2, RS3, RS4,
     *                  SAFETY, SQRRMC, STBRAD, T, TANANG, TINY, TOLD,
     *                  TOOSML
      INTEGER           FLSTP, GNFCN, LSTSTG, MAXTRY, NFCN, NSEC, OKSTP,
     *                  ORDER, OUTCH, PRZERR, PRZERS, PRZSTG, PRZY,
     *                  PRZYNU, PRZYP, SVNFCN
      LOGICAL           ERASFL, ERASON, FIRST, FSAL, LAST
C     .. Local Scalars ..
      DOUBLE PRECISION  DIFF, ERRMAX, HMIN, HSEC, MXERLC, TSEC, ZLERR,
     *                  ZTEST1, ZTEST2
      INTEGER           ISTEP, L, LEVEL
      LOGICAL           LDUMMY, MAIN
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMARR(1)
C     .. External Subroutines ..
      EXTERNAL          D02PDZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX
C     .. Common blocks ..
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /ED02PD/TOOSML, COST, SAFETY, EXPON, STBRAD,
     *                  TANANG, RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG,
     *                  MAXTRY, NSEC, FSAL
      COMMON            /FD02PD/MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY,
     *                  PRZYP, PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
      COMMON            /GD02PD/MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC,
     *                  TINY, OUTCH
C     .. Save statement ..
      SAVE              /BD02PD/, /ED02PD/, /FD02PD/, /GD02PD/
C     .. Executable Statements ..
      TSEC = T - HOLD
      HSEC = HOLD/DBLE(NSEC)
      HMIN = MAX(TINY,TOOSML*MAX(ABS(TSEC),ABS(T)))
      IF (ABS(HSEC).LT.HMIN) THEN
         IER = 6
         GO TO 120
      END IF
      ZTEST1 = TOL/DBLE(NSEC)
      ZTEST2 = TOL/TEN
      LEVEL = 0
C
C  The subroutine D02PDZ is used to take a step. In its use in the
C  primary integration provision is made for getting on scale in the
C  first step. In this situation only the subroutine might reduce the
C  step size. By setting MAIN = .FALSE., the subroutine will take a
C  step of the size input. In this use of the subroutine, all items of
C  the call list appearing after MAIN are dummy variables.
C
C  Perform secondary integration.
C
      MAIN = .FALSE.
      LDUMMY = .FALSE.
      DO 60 ISTEP = 1, NSEC
C
C  Take a step.
C
         CALL D02PDZ(F,NEQ,TSEC,ZY,ZYP,ZSTAGE,ZTEST1,HSEC,WEIGHT,ZYNEW,
     *               ZERRES,ZLERR,MAIN,DUMMY,DUMARR,LDUMMY)
C
C  The primary integration is using a step size of HUSED and the
C  secondary integration is using the smaller step size
C  HSEC = HUSED/NSEC.  If steps of this size were taken from the same
C  starting point and the asymptotic behavior were evident, the smaller
C  step size would result in a local error that is considerably smaller,
C  namely by a factor of 1/(NSEC**(ORDER+1)). If the two approximate
C  solutions are close and TOLR is neither too large nor too small, this
C  should be approximately true.  The step size is chosen in the primary
C  integration so that the local error ERR is no larger than TOLR. The
C  local error, ZLERR, of the secondary integration is compared to TOLR
C  in an attempt to diagnose a secondary integration that is not rather
C  more accurate than the primary integration.
C
         IF (ZLERR.GE.ZTEST1) THEN
            LEVEL = 2
         ELSE IF (ZLERR.GT.ZTEST2) THEN
            LEVEL = LEVEL + 1
         END IF
         IF (LEVEL.GE.2) THEN
            IER = 6
            GO TO 120
         END IF
C
C  Advance TSEC and the dependent variables ZY(*) and ZYP(*).
         TSEC = T - DBLE(NSEC-ISTEP)*HSEC
         DO 20 L = 1, NEQ
            ZY(L) = ZYNEW(L)
   20    CONTINUE
C
         IF (FSAL) THEN
C
C  When FSAL = .TRUE., the derivative ZYP(*) is the last stage of the
C  step.
C
            DO 40 L = 1, NEQ
               ZYP(L) = ZSTAGE(L,LSTSTG)
   40       CONTINUE
         ELSE
C
C  Call F to evaluate ZYP(*).
            CALL F(TSEC,ZY,ZYP)
            GNFCN = GNFCN + 1
         END IF
C
   60 CONTINUE
C
C  Update the maximum error seen, MAXERR, and its location, LOCMAX.
C  Use local variables ERRMAX and MXERLC.
C
      ERRMAX = MAXERR
      MXERLC = LOCMAX
      DO 80 L = 1, NEQ
         DIFF = ABS(ZY(L)-Y(L))/WEIGHT(L)
         IF (DIFF.GT.ERRMAX) THEN
            ERRMAX = DIFF
            MXERLC = T
         END IF
   80 CONTINUE
C
C  If the global error is greater than 0.1D0, the solutions have
C  diverged so far that comparing them may not provide a reliable
C  estimate of the global error. The test is made before ZERROR(*)
C  and MAXERR, LCMXER are updated so that on a failure, they refer
C  to the last reliable results.
C
      IF (ERRMAX.GT.PT1) THEN
         IER = 6
         GO TO 120
      ELSE
         MAXERR = ERRMAX
         LOCMAX = MXERLC
         DO 100 L = 1, NEQ
            DIFF = ABS(ZY(L)-Y(L))/WEIGHT(L)
            ZERROR(L) = ZERROR(L) + DIFF**2
  100    CONTINUE
         IER = 0
      END IF
C
C  Exit point for D02PDW
  120 CONTINUE
C
      RETURN
      END
