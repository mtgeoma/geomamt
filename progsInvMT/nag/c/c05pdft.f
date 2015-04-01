      SUBROUTINE C05PDF(IREVCM,N,X,FVEC,FJAC,LDFJAC,XTOL,DIAG,MODE,
     *                  FACTOR,R,LR,QTF,W,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     **********
C
C     SUBROUTINE C05PDF
C
C     The purpose of C05PDF is to interface to C05PDZ.
C     The latter is based upon MINPACK routine HYBRJ.
C
C     **********
C
C     Revised to output explanatory messages.
C     P.J.D. Mayes, NAG Central Office, December 1987.
C     Revised to reverse communication.
C     M.S. Derakhshan, NAG Central Office, September 1988.
C
C     ..Subroutine Arguments ..
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05PDF')
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FACTOR, XTOL
      INTEGER           IFAIL, IREVCM, LDFJAC, LR, MODE, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DIAG(N), FJAC(LDFJAC,N), FVEC(N), QTF(N), R(LR),
     *                  W(N,4), X(N)
C     .. Local Scalars ..
      INTEGER           INFO, LR1, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05PDZ
C     .. Save statement ..
      SAVE
C     .. Executable Statements ..
C
      IF (IREVCM.NE.0) GO TO 20
      INFO = 1
      LR1 = N*(N+1)/2
      NREC = 0
      IF (N.LE.0) THEN
         WRITE (P01REC,FMT=99999) N
         NREC = 1
         GO TO 40
      ELSE IF (LDFJAC.LT.N) THEN
         WRITE (P01REC,FMT=99998) LDFJAC, N
         NREC = 2
         GO TO 40
      ELSE IF (XTOL.LT.ZERO) THEN
         WRITE (P01REC,FMT=99997) XTOL
         NREC = 1
         GO TO 40
      ELSE IF (FACTOR.LE.ZERO) THEN
         WRITE (P01REC,FMT=99996) FACTOR
         NREC = 1
         GO TO 40
      ELSE IF (LR.LT.LR1) THEN
         WRITE (P01REC,FMT=99995) LR, LR1
         NREC = 2
         GO TO 40
      END IF
   20 IF (IREVCM.LT.0 .OR. IREVCM.GT.3) THEN
         INFO = 2
         WRITE (P01REC,FMT=99993) IREVCM
         NREC = 1
         GO TO 40
      END IF
      CONTINUE
      CALL C05PDZ(IREVCM,N,X,FVEC,FJAC,LDFJAC,XTOL,DIAG,MODE,FACTOR,
     *            INFO,R,LR,QTF,W(1,1),W(1,2),W(1,3),W(1,4))
      IF (IREVCM.GT.0) RETURN
      IF (INFO.EQ.1) THEN
         WRITE (P01REC,FMT=99992)
         NREC = 1
      ELSE IF (INFO.EQ.3) THEN
         WRITE (P01REC,FMT=99994) XTOL
         NREC = 2
      ELSE IF (INFO.EQ.4) THEN
         WRITE (P01REC,FMT=99991)
         NREC = 2
      ELSE IF (INFO.EQ.5) THEN
         WRITE (P01REC,FMT=99990)
         NREC = 2
      END IF
   40 IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N must be greater than 0: N = ',I16)
99998 FORMAT (' ** On entry, LDFJAC must be at least N:',/' ** LDFJAC ',
     *       '= ',I16,',    N = ',I16)
99997 FORMAT (' ** On entry, XTOL must be at least 0.0: XTOL = ',1P,
     *       D13.5)
99996 FORMAT (' ** On entry, FACTOR must be greater than 0.0: FACTOR = '
     *       ,1P,D13.5)
99995 FORMAT (' ** On entry, LR must be at least N*(N+1)/2:',/' ** LR ',
     *       '= ',I16,',    N*(N+1)/2 = ',I16)
99994 FORMAT (' ** No further improvement in the solution is possible.',
     *       /' ** XTOL is too small: XTOL = ',1P,D13.5)
99993 FORMAT (' ** On entry, IREVCM was not set to 0,1,2,3: IREVCM = ',
     *       I16)
99992 FORMAT (' ** On entry, MODE=2 and DIAG contained a non-positive ',
     *       'element.')
99991 FORMAT (' ** The iteration is not making good progress as measur',
     *       'ed',/' ** by the improvement from the last 5 Jacobian ev',
     *       'aluations')
99990 FORMAT (' ** The iteration is not making good progress as measur',
     *       'ed',/' ** by the improvement from the last 10 iterations '
     *       )
      END
