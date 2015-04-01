      SUBROUTINE C05NBF(FCN,N,X,FVEC,TOL,WA,LWA,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05NBF (based on MINPACK routine HYBRD1)
C
C     The purpose of C05NBF is to find a zero of a system of
C     N nonlinear functions in N variables by a modification
C     of the Powell Hybrid method. This is done by using the
C     more general nonlinear equation solver C05NCF. The user
C     must provide a subroutine which calculates the functions.
C     The Jacobian is then calculated by a forward-difference
C     approximation.
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C     **********
C
C     Revised to output explanatory messages.
C     P.J.D. Mayes, NAG Central Office, December 1987
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05NBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, LWA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FVEC(N), WA(LWA), X(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSFCN, FACTOR, ONE, XTOL, ZERO
      INTEGER           INDEX, INFO, J, LR, LR1, MAXFEV, ML, MODE, MU,
     *                  NFEV, NPRINT, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05NCS
C     .. Data statements ..
      DATA              FACTOR, ONE, ZERO/1.0D2, 1.0D0, 0.0D0/
C     .. Executable Statements ..
C
C     Check the input parameters for errors.
C
      INFO = 1
      LR1 = N*(3*N+13)/2
      NREC = 0
      IF (N.LE.0) THEN
         WRITE (P01REC,FMT=99999) N
         NREC = 1
         GO TO 40
      ELSE IF (TOL.LT.ZERO) THEN
         WRITE (P01REC,FMT=99998) TOL
         NREC = 1
         GO TO 40
      ELSE IF (LWA.LT.LR1) THEN
         WRITE (P01REC,FMT=99997) LWA, LR1
         NREC = 2
         GO TO 40
      END IF
C
C     Call C05NCS.
C
      MAXFEV = 200*(N+1)
      XTOL = TOL
      ML = N - 1
      MU = N - 1
      EPSFCN = ZERO
      MODE = 2
      DO 20 J = 1, N
         WA(J) = ONE
   20 CONTINUE
      NPRINT = 0
      LR = (N*(N+1))/2
      INDEX = 6*N + LR
      CALL C05NCS(FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,WA(1),MODE,
     *            FACTOR,NPRINT,INFO,NFEV,WA(INDEX+1),N,WA(6*N+1),LR,
     *            WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))
C
C     Check for errors on exit from C05NCS
C
      IF (INFO.LT.0) THEN
         P01REC(1) = ' ** User set IFLAG negative in FCN'
         NREC = 1
      ELSE IF (INFO.EQ.2) THEN
         P01REC(1) =
     *             ' ** There have been at least 200*(N+1) calls to FCN'
         NREC = 1
      ELSE IF (INFO.EQ.3) THEN
         WRITE (P01REC,FMT=99996) XTOL
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
99998 FORMAT (' ** On entry, XTOL must be at least 0.0: XTOL = ',1P,
     *       D13.5)
99997 FORMAT (' ** On entry, LWA must be at least N*(3*N+13)/2:',/' **',
     *       ' LWA = ',I16,',    N*(3*N+13)/2 = ',I16)
99996 FORMAT (' ** No further improvement in the solution is possible.',
     *       /' ** XTOL is too small: XTOL = ',1P,D13.5)
      END
