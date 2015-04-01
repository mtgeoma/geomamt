      SUBROUTINE C05NDF(IREVCM,N,X,FVEC,XTOL,ML,MU,EPSFCN,DIAG,MODE,
     *                  FACTOR,FJAC,LDFJAC,R,LR,QTF,W,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     **********
C
C     SUBROUTINE C05NDF
C
C     The purpose of C05NDF is to interface to C05NDS.
C     The latter is based on MINPACK routine HYBRD.
C
C     **********
C
C     Revised to output explanatory messages.
C     P.J.D. Mayes, NAG Central Office, December 1987.
C     Revised to reverse communication.
C     M.S. Derakhshan, NAG Central Office, September 1988.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05NDF')
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSFCN, FACTOR, XTOL
      INTEGER           IFAIL, IREVCM, LDFJAC, LR, ML, MODE, MU, N
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
      EXTERNAL          C05NDS
C     .. Save statement ..
      SAVE
C     .. Executable Statements ..
      IF (IREVCM.NE.0) GO TO 20
      INFO = 1
      LR1 = N*(N+1)/2
      NREC = 0
      IF (N.LE.0) THEN
         WRITE (P01REC,FMT=99997) N
         NREC = 1
         GO TO 40
      ELSE IF (LDFJAC.LT.N) THEN
         WRITE (P01REC,FMT=99996) LDFJAC, N
         NREC = 2
         GO TO 40
      ELSE IF (ML.LT.0) THEN
         WRITE (P01REC,FMT=99999) ML
         NREC = 1
         GO TO 40
      ELSE IF (MU.LT.0) THEN
         WRITE (P01REC,FMT=99998) MU
         NREC = 1
         GO TO 40
      ELSE IF (XTOL.LT.ZERO) THEN
         WRITE (P01REC,FMT=99995) XTOL
         NREC = 1
         GO TO 40
      ELSE IF (FACTOR.LE.ZERO) THEN
         WRITE (P01REC,FMT=99994) FACTOR
         NREC = 1
         GO TO 40
      ELSE IF (LR.LT.LR1) THEN
         WRITE (P01REC,FMT=99993) LR, LR1
         NREC = 2
         GO TO 40
      END IF
   20 IF (IREVCM.LT.0 .OR. IREVCM.GT.2) THEN
         INFO = 2
         WRITE (P01REC,FMT=99991) IREVCM
         NREC = 1
         GO TO 40
      END IF
      CONTINUE
      CALL C05NDS(IREVCM,N,X,FVEC,XTOL,ML,MU,EPSFCN,DIAG,MODE,FACTOR,
     *            INFO,FJAC,LDFJAC,R,LR,QTF,W(1,1),W(1,2),W(1,3),W(1,4))
      IF (IREVCM.GT.0) RETURN
      IF (INFO.EQ.1) THEN
         WRITE (P01REC,FMT=99990)
         NREC = 1
      ELSE IF (INFO.EQ.3) THEN
         WRITE (P01REC,FMT=99992) XTOL
         NREC = 2
      ELSE IF (INFO.EQ.4) THEN
         WRITE (P01REC,FMT=99989)
         NREC = 2
      ELSE IF (INFO.EQ.5) THEN
         WRITE (P01REC,FMT=99988)
         NREC = 2
      END IF
   40 IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, ML must be at least 0: ML = ',I16)
99998 FORMAT (' ** On entry, MU must be at least 0: MU = ',I16)
99997 FORMAT (' ** On entry, N must be greater than 0: N = ',I16)
99996 FORMAT (' ** On entry, LDFJAC must be at least N:',/' ** LDFJAC ',
     *       '= ',I16,',    N = ',I16)
99995 FORMAT (' ** On entry, XTOL must be at least 0.0: XTOL = ',1P,
     *       D13.5)
99994 FORMAT (' ** On entry, FACTOR must be greater than 0.0: FACTOR = '
     *       ,1P,D13.5)
99993 FORMAT (' ** On entry, LR must be at least N*(N+1)/2:',/' ** LR ',
     *       '= ',I16,',    N*(N+1)/2 = ',I16)
99992 FORMAT (' ** No further improvement in the solution is possible.',
     *       /' ** XTOL is too small: XTOL = ',1P,D13.5)
99991 FORMAT (' ** On entry, IREVCM was not set to 0,1,2: IREVCM = ',
     *       I16)
99990 FORMAT (' ** On entry, MODE=2 and DIAG contained a non-positive ',
     *       'element.')
99989 FORMAT (' ** The iteration is not making good progress as measur',
     *       'ed',/' ** by the improvement from the last 5 Jacobian ev',
     *       'aluations')
99988 FORMAT (' ** The iteration is not making good progress as measur',
     *       'ed',/' ** by the improvement from the last 10 iterations')
      END
