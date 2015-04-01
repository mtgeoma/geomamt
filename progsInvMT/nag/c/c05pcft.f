      SUBROUTINE C05PCF(FCN,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,DIAG,MODE,
     *                  FACTOR,NPRINT,NFEV,NJEV,R,LR,QTF,W,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05PCF
C
C     The purpose of C05PCF is to interface to C05PCZ.
C     The latter is based upon MINPACK routine HYBRJ.
C
C     **********
C
C     Revised to output explanatory messages.
C     P.J.D. Mayes, NAG Central Office, December 1987.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05PCF')
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FACTOR, XTOL
      INTEGER           IFAIL, LDFJAC, LR, MAXFEV, MODE, N, NFEV, NJEV,
     *                  NPRINT
C     .. Array Arguments ..
      DOUBLE PRECISION  DIAG(N), FJAC(LDFJAC,N), FVEC(N), QTF(N), R(LR),
     *                  W(N,4), X(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      INTEGER           INFO, LR1, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05PCZ
C     .. Executable Statements ..
      INFO = 1
      LR1 = N*(N+1)/2
      NREC = 0
      IF (N.LE.0) THEN
         WRITE (P01REC,FMT=99999) N
         NREC = 1
         GO TO 20
      ELSE IF (LDFJAC.LT.N) THEN
         WRITE (P01REC,FMT=99998) LDFJAC, N
         NREC = 2
         GO TO 20
      ELSE IF (XTOL.LT.ZERO) THEN
         WRITE (P01REC,FMT=99997) XTOL
         NREC = 1
         GO TO 20
      ELSE IF (MAXFEV.LE.0) THEN
         WRITE (P01REC,FMT=99996) MAXFEV
         NREC = 1
         GO TO 20
      ELSE IF (FACTOR.LE.ZERO) THEN
         WRITE (P01REC,FMT=99995) FACTOR
         NREC = 1
         GO TO 20
      ELSE IF (LR.LT.LR1) THEN
         WRITE (P01REC,FMT=99994) LR, LR1
         NREC = 2
         GO TO 20
      END IF
      CALL C05PCZ(FCN,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,DIAG,MODE,FACTOR,
     *            NPRINT,INFO,NFEV,NJEV,R,LR,QTF,W(1,1),W(1,2),W(1,3),
     *            W(1,4))
      IF (INFO.LT.0) THEN
         P01REC(1) = ' ** User set IFLAG negative in FCN'
         NREC = 1
      ELSE IF (INFO.EQ.1) THEN
         P01REC(1) =
     * ' ** On entry, MODE=2 and DIAG contained a non-positive element.'
         NREC = 1
      ELSE IF (INFO.EQ.2) THEN
         WRITE (P01REC,FMT=99993) MAXFEV
         NREC = 1
      ELSE IF (INFO.EQ.3) THEN
         WRITE (P01REC,FMT=99992) XTOL
         NREC = 2
      ELSE IF (INFO.EQ.4) THEN
         P01REC(1) =
     *       ' ** The iteration is not making good progress as measured'
         P01REC(2) =
     *     ' ** by the improvement from the last 5 Jacobian evaluations'
         NREC = 2
      ELSE IF (INFO.EQ.5) THEN
         P01REC(1) =
     *       ' ** The iteration is not making good progress as measured'
         P01REC(2) =
     *              ' ** by the improvement from the last 10 iterations'
         NREC = 2
      END IF
   20 IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N must be greater than 0: N = ',I16)
99998 FORMAT (' ** On entry, LDFJAC must be at least N:',/' ** LDFJAC ',
     *       '= ',I16,',    N = ',I16)
99997 FORMAT (' ** On entry, XTOL must be at least 0.0: XTOL = ',1P,
     *       D13.5)
99996 FORMAT (' ** On entry, MAXFEV must be greater than 0: MAXFEV = ',
     *       I16)
99995 FORMAT (' ** On entry, FACTOR must be greater than 0.0: FACTOR = '
     *       ,1P,D13.5)
99994 FORMAT (' ** On entry, LR must be at least N*(N+1)/2:',/' ** LR ',
     *       '= ',I16,',    N*(N+1)/2 = ',I16)
99993 FORMAT (' ** There have been at least MAXFEV calls to FCN: MAXFE',
     *       'V = ',I16)
99992 FORMAT (' ** No further improvement in the solution is possible.',
     *       /' ** XTOL is too small: XTOL = ',1P,D13.5)
      END
