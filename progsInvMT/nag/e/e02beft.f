      SUBROUTINE E02BEF(START,M,X,Y,W,S,NEST,N,K,C,FP,WRK,LWRK,IWRK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02BEF')
      INTEGER           K0, MAXIT
      DOUBLE PRECISION  TOL, ZERO
      PARAMETER         (K0=3,MAXIT=20,TOL=1.0D-3,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FP, S
      INTEGER           IFAIL, LWRK, M, N, NEST
      CHARACTER*1       START
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NEST), K(NEST), W(M), WRK(LWRK), X(M), Y(M)
      INTEGER           IWRK(NEST)
C     .. Local Scalars ..
      DOUBLE PRECISION  TS, XB, XE
      INTEGER           I, IA, IB, IER, IFP, IG, IOPT, IQ, IZ, K1, K2,
     *                  LWEST, NMIN, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02BEZ
C     .. Executable Statements ..
C     Set up the boundaries of the approximation interval.
      XB = X(1)
      XE = X(M)
C     Before starting computations a data check is made. If the
C     input data are invalid, control is immediately repassed
C     to the calling program.
      IER = 1
      K1 = K0 + 1
      K2 = K1 + 1
      NMIN = 2*K1
      LWEST = (M+K1)*K1 + NEST*(7+3*K0) + K2*K2
      NREC = 1
      IF (S.GT.ZERO .AND. S.LT.X02AJF()) THEN
C        If S is smaller than machine precision, assume it is zero.
         TS = ZERO
      ELSE
         TS = S
      END IF
      IF (START.EQ.'C' .OR. START.EQ.'c') THEN
         IOPT = 0
      ELSE IF (START.EQ.'W' .OR. START.EQ.'w') THEN
         IOPT = 1
      ELSE
         IOPT = -1
      END IF
      IF (IOPT.EQ.-1) THEN
         WRITE (REC,FMT=99999) START
      ELSE IF (M.LT.K1) THEN
         WRITE (REC,FMT=99998) M
      ELSE IF (TS.LT.ZERO) THEN
         WRITE (REC,FMT=99997) TS
      ELSE IF (TS.EQ.ZERO .AND. NEST.LT.M+4) THEN
         WRITE (REC,FMT=99996) M + 4, NEST
      ELSE IF (NEST.LT.NMIN) THEN
         WRITE (REC,FMT=99995) NEST
      ELSE IF (LWRK.LT.LWEST) THEN
         WRITE (REC,FMT=99994) LWEST, LWRK
      ELSE IF (W(1).LE.ZERO) THEN
         WRITE (REC,FMT=99993) 1
      ELSE
         DO 20 I = 2, M
            IF (W(I).LE.ZERO) THEN
               IER = 2
               NREC = 2
               WRITE (REC,FMT=99993) I, W(I)
               GO TO 40
            ELSE IF (X(I-1).GE.X(I)) THEN
               IER = 3
               NREC = 2
               WRITE (REC,FMT=99992) I, X(I-1), X(I)
               GO TO 40
            END IF
   20    CONTINUE
C        Partition the working space and determine the spline
C        approximation.
         IER = 0
         NREC = 0
         IFP = 1
         IZ = IFP + NEST
         IA = IZ + NEST
         IB = IA + (NEST+K1)*K1
         IG = IB + NEST*K2
         IQ = IG + (NEST+K2)*K2
         CALL E02BEZ(IOPT,X,Y,W,M,XB,XE,K0,TS,NEST,TOL,MAXIT,K1,K2,N,K,
     *               C,FP,WRK(IFP),WRK(IZ),WRK(IA),WRK(IB),WRK(IG),
     *               WRK(IQ),IWRK,IER)
         IF (IER.EQ.4) THEN
            NREC = 2
            WRITE (REC,FMT=99991) NEST, TS
         ELSE IF (IER.EQ.5) THEN
            NREC = 2
            WRITE (REC,FMT=99990) TS
         END IF
      END IF
   40 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, START .ne. ''W'' or ''C'': START = ''',A,
     *       '''')
99998 FORMAT (1X,'** On entry, M .lt. 4: M =',I16)
99997 FORMAT (1X,'** On entry, S .lt. zero: S =',1P,D13.5)
99996 FORMAT (1X,'** On entry, NEST .lt.',I16,' when S = zero: NEST =',
     *       I16)
99995 FORMAT (1X,'** On entry, NEST .lt. 8: NEST =',I16)
99994 FORMAT (1X,'** On entry, LWRK .lt.',I16,': LWRK =',I16)
99993 FORMAT (1X,'** On entry, W(R) .le. zero for R =',I16,':',/4X,'W(',
     *       'R) =',1P,D13.5,' .')
99992 FORMAT (1X,'** On entry, X(R-1) .ge. X(R) for R =',I16,':',/4X,
     *       'X(R-1), X(R) =',1P,2D13.5,' .')
99991 FORMAT (1X,'** The number of knots needed is greater than NEST: ',
     *       'NEST =',I8,' .',/4X,'Possibly S is too small: S =',1P,
     *       D13.5,' .')
99990 FORMAT (1X,'** The iterative process has failed to converge.',/4X,
     *       'Possibly S is too small: S =',1P,D13.5,' .')
      END
