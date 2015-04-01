      SUBROUTINE C05PBF(FCN,N,X,FVEC,FJAC,LDFJAC,TOL,WA,LWA,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05PBF (based on MINPACK routine HYBRJ1)
C
C     The purpose of C05PBF is to find a zero of a system of
C     N nonlinear functions in N variables by a modification
C     of the Powell Hybrid method. This is done by using the
C     more general nonlinear equation solver C05PCF. The user
C     must provide a subroutine which calculates the functions
C     and the Jacobian.
C
C     **********
C
C     Revised to output explanatory messages.
C     P.J.D. Mayes, NAG Central Office, December 1987.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05PBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, LDFJAC, LWA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LDFJAC,N), FVEC(N), WA(LWA), X(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  FACTOR, ONE, XTOL, ZERO
      INTEGER           INFO, J, LR, LR1, MAXFEV, MODE, NFEV, NJEV,
     *                  NPRINT, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05PCZ
C     .. Data statements ..
      DATA              FACTOR, ONE, ZERO/1.0D2, 1.0D0, 0.0D0/
C     .. Executable Statements ..
      INFO = 1
C
C     Check the input parameters for errors.
C
      INFO = 1
      LR1 = N*(N+13)/2
      NREC = 0
      IF (N.LE.0) THEN
         WRITE (P01REC,FMT=99999) N
         NREC = 1
         GO TO 40
      ELSE IF (LDFJAC.LT.N) THEN
         WRITE (P01REC,FMT=99998) LDFJAC, N
         NREC = 2
         GO TO 40
      ELSE IF (TOL.LT.ZERO) THEN
         WRITE (P01REC,FMT=99997) TOL
         NREC = 1
         GO TO 40
      ELSE IF (LWA.LT.LR1) THEN
         WRITE (P01REC,FMT=99996) LWA, LR1
         NREC = 2
         GO TO 40
      END IF
C
C     Call C05PCZ.
C
      MAXFEV = 100*(N+1)
      XTOL = TOL
      MODE = 2
      DO 20 J = 1, N
         WA(J) = ONE
   20 CONTINUE
      NPRINT = 0
      LR = (N*(N+1))/2
      CALL C05PCZ(FCN,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,WA(1),MODE,
     *            FACTOR,NPRINT,INFO,NFEV,NJEV,WA(6*N+1),LR,WA(N+1),
     *            WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))
      IF (INFO.LT.0) THEN
         P01REC(1) = ' ** User set IFLAG negative in FCN'
         NREC = 1
      ELSE IF (INFO.EQ.2) THEN
         P01REC(1) =
     *             ' ** There have been at least 100*(N+1) calls to FCN'
         NREC = 1
      ELSE IF (INFO.EQ.3) THEN
         WRITE (P01REC,FMT=99995) XTOL
         NREC = 2
      ELSE IF (INFO.EQ.4 .OR. INFO.EQ.5) THEN
         P01REC(1) = ' ** The iteration is not making good progress'
         NREC = 1
         INFO = 4
      END IF
   40 IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N must be greater than 0: N = ',I16)
99998 FORMAT (' ** On entry, LDFJAC must be at least N:',/' ** LDFJAC ',
     *       '= ',I16,', N = ',I16)
99997 FORMAT (' ** On entry, XTOL must be at least 0.0: XTOL = ',1P,
     *       D13.5)
99996 FORMAT (' ** On entry, LWA must be at least N*(N+13)/2:',/' ** L',
     *       'WA = ',I16,',    N*(N+13)/2 = ',I16)
99995 FORMAT (' ** No further improvement in the solution is possible.',
     *       /' ** XTOL is too small: XTOL = ',1P,D13.5)
      END
