      SUBROUTINE S14BAF(A,X,TOL,P,Q,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF, QUARTR, THRTWO, TWO, THREE,
     *                  FOUR, EIGHT, TWENTY, PSEVEN, ONEP4
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,QUARTR=0.25D0,
     *                  THRTWO=1.5D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0,
     *                  EIGHT=8.0D0,TWENTY=20.0D0,PSEVEN=0.7D0,
     *                  ONEP4=1.4D0)
      INTEGER           MAXIT, NTERMS
      PARAMETER         (MAXIT=600,NTERMS=22)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='S14BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, P, Q, TOL, X
      INTEGER           IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION  ALG, ALGP1, ALGS, ALPHA, ALX, AP1, BOT, EPS,
     *                  EPS1, GA, PP, QQ, RHO, RR, SAFE, SS, SUM, TERM,
     *                  TT, U, UNDFL, V, XMA, XPA, Y
      INTEGER           IERR, K, NREC, TIFAIL
C     .. Local Arrays ..
      DOUBLE PRECISION  C(NTERMS)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S14AAF, S14ABF, S14BAZ, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          S14AAF, S14ABF, S14BAZ, X02AJF, X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, SQRT
C     .. Data statements ..
      DATA              C/5.7721566490153286060651209008D-1,
     *                  -6.5587807152025388107701951515D-1,
     *                  -4.200263503409523552900393488D-2,
     *                  1.665386113822914895017007951D-1,
     *                  -4.21977345555443367482083013D-2,
     *                  -9.6219715278769735621149217D-3,
     *                  7.21894324666309954230950103D-3,
     *                  -1.16516759185906511211971D-3,
     *                  -2.15241674114950972815730D-4,
     *                  1.2805028238811618615320D-4,
     *                  -2.013485478078823865569D-5,
     *                  -1.2504934821426706573D-6,
     *                  1.1330272319816958824D-6,
     *                  -2.056338416977607103D-7, 6.1160951044814158D-9,
     *                  5.0020076444692229D-9, -1.181274570487020D-9,
     *                  1.04342671169110D-10, 7.782263439905D-12,
     *                  -3.696805618642D-12, 5.1003702875D-13,
     *                  -2.058326054D-14/
C     .. Executable Statements ..
C
C     Let  GAMMA(A)  denote the gamma function and  GAM(A,X)  the
C     (complementary) incomplete gamma function,
C
C     GAM(A,X)=integral from T=X to T=infinity of EXP(-T)*T**(A-1).
C
C     Let  GAMSTAR(A,X)  denote Tricomi's form of the incomplete gamma
C     function, which for A.gt.0. is defined by
C
C     GAMSTAR(A,X)=(X**(-A)/GAMMA(A))*integral from T=0 to T=X of
C                EXP(-T)*T**(A-1).
C
C     For the purpose of this subroutine, these functions are normalized
C     as follows:
C
C     P(A,X)  =   (X**A)*GAMSTAR(A,X).
C
C     Q(A,X)  =   GAM(A,X)/GAMMA(A),
C
C     The program below attempts to evaluate  P(A,X)  and  Q(A,X),
C     both to an accuracy of TOL significant decimal digits, for
C     positive A and nonnegative X. There are (rare) instances in
C     which the accuracy attained is somewhat less than the accuracy
C     specified. The discrepancy, however, should never exceed one or
C     two (decimal) orders of accuracy.
C     The functions are evaluated by means of the Maclaurin expansion,
C     Taylor expansion, or Legendre's continued fraction, depending on
C     the values of A and X. However, when A .ge. 20.0 and
C     0.7*A .le. X .le. 1.4*A, the asymptotic expansion of Temme is
C     used for greater efficiency.
C
C     Arguments:
C         A - The first argument of P and Q.
C         X - The second argument of P and Q.
C       TOL - The number of correct significant decimal digits
C             desired in the results.
C         P - An output variable returning the value of P(A,X).
C         Q - An output variable returning the value of Q(A,X).
C     IFAIL - Error flag. The values of IFAIL have these meanings:
C            0 - Successful exit.
C            1 - Illegal negative or zero argument A. The routine
C                exits with the value zero for P and Q.
C            2 - Illegal negative argument X. The routine exits
C                with the value zero for P and Q.
C            3 - Convergence fails within MAXIT (=600) iterations,
C                either in Taylor's series or in Legendre's continued
C                fraction. Reason unknown. The computation is
C                aborted and the routine exits with the value zero
C                for P and Q.
C
C     The data declaration contains the successive coefficients in the
C     Maclaurin expansion of (1/GAMMA(A+1))-1.
C     Values of these coefficients (to 31 decimal places) can be found
C     in Table 5 of J.W.Wrench,jr.:
C     Concerning Two Series for the Gamma Function  , Math. Comput.
C     22, 1968, 617-626.
C
C     This routine is derived from ACM Algorithm 542.
C
C     References:
C     W. Gautschi,   'A Computational Procedure for Incomplete Gamma
C     Functions', ACM Trans. Math. Software, Vol. 5, No. 4, 1979,
C     466-481.
C     N. M. Temme,   'On the computation of the incomplete gamma
C     functions for large values of the parameters', proceedings of
C     Shrivenham conference 'Algorithms for Approximation',
C     ed. J.C.Mason and M.G.Cox, Oxford University Press 1987.
C
      P = ZERO
      Q = ZERO
      NREC = 0
      IERR = 0
      IF (A.LE.ZERO) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) A
      ELSE IF (X.LT.ZERO) THEN
         IERR = 2
         NREC = 1
         WRITE (REC,FMT=99998) X
      ELSE IF (X.EQ.ZERO) THEN
         P = ZERO
         Q = ONE
      ELSE
C        1.0E-18 is the accuracy of constants in this routine and
C        auxiliaries.
         EPS = MAX(X02AJF(),1.0D-18)
         IF (TOL.GT.EPS .AND. TOL.LE.ONE) EPS = TOL
         SAFE = X02AMF()
         IF (A.GE.TWENTY .AND. PSEVEN*A.LE.X .AND. X/ONEP4.LE.A) THEN
C           Use the asymptotic expansion of Temme.
            UNDFL = LOG(SAFE)
            IF (A.LE.X) THEN
               Q = S14BAZ(A,X,.FALSE.,EPS,UNDFL)
               P = ONE - Q
            ELSE
               P = S14BAZ(A,X,.TRUE.,EPS,UNDFL)
               Q = ONE - P
            END IF
         ELSE
            ALX = LOG(X)
            IF (X.LT.QUARTR) THEN
               ALPHA = LOG(HALF)/ALX
            ELSE
               ALPHA = X + QUARTR
            END IF
            BOT = LOG(SAFE)
            EPS1 = EPS/100
            AP1 = A + ONE
C
C           Evaluation of the logarithm of GAMMA(A+1.0).
C
            TIFAIL = 1
            ALGP1 = S14ABF(AP1,TIFAIL)
            IF (TIFAIL.EQ.2) THEN
C              ln gamma(a+1) overflows. P and Q are 1.0 and 0.0,
C              or vice versa, to machine precision.
               IF (A.GT.X) THEN
                  P = ZERO
                  Q = ONE
               ELSE
                  P = ONE
                  Q = ZERO
               END IF
            ELSE IF (A.GT.ALPHA) THEN
C
C              Evaluation of P(A,X) for A.gt.ALPHA(X) by Taylor
C              expansion.
C
               TERM = ONE
               SUM = ONE
               K = 0
   20          K = K + 1
               IF (K.GT.MAXIT) GO TO 120
               TERM = X*TERM/(A+K)
               SUM = SUM + TERM
               IF (ABS(TERM).GT.EPS*SUM) GO TO 20
               ALGS = A*ALX - X + LOG(SUM) - ALGP1
               IF (ALGS.LE.BOT) THEN
                  P = ZERO
               ELSE
                  P = EXP(ALGS)
               END IF
               Q = ONE - P
            ELSE IF (X.GT.THRTWO) THEN
C
C              Evaluation of Q(A,X) for X.gt.1.5 and A.le.ALPHA(X) by
C              means of the Legendre continued fraction.
C
               XPA = X + ONE - A
               XMA = X - ONE - A
               IF (XPA.GT.ONE/SQRT(SAFE)) THEN
C                 XPA*XMA would overflow, but P = 1.0 to machine
C                 precision.
                  P = ONE
                  Q = ZERO
               ELSE
                  PP = ZERO
                  QQ = XPA*XMA
                  RR = FOUR*XPA
                  SS = -A + ONE
                  TERM = ONE
                  SUM = ONE
                  RHO = ZERO
                  K = 1
   40             K = K + 1
                  IF (K.GT.MAXIT) GO TO 120
                  PP = PP + SS
                  QQ = QQ + RR
                  RR = RR + EIGHT
                  SS = SS + TWO
                  TT = PP*(ONE+RHO)
                  RHO = TT/(QQ-TT)
                  TERM = RHO*TERM
                  SUM = SUM + TERM
                  IF (ABS(TERM).GT.EPS*SUM) GO TO 40
                  ALG = A*ALX - X + LOG(A*SUM/XPA) - ALGP1
                  IF (ALG.LE.BOT) THEN
                     Q = ZERO
                  ELSE
                     Q = EXP(ALG)
                  END IF
                  P = ONE - Q
               END IF
            ELSE
C
C              Direct evaluation of Q(A,X) and P(A,X) for X.le.1.5
C              and A.le.ALPHA(X).
C
               IF (A.GT.HALF) THEN
                  TIFAIL = -1
                  U = S14AAF(A,TIFAIL) - (X**A)/A
               ELSE
C                 NTERMS = 22 is sufficient for approximately 18
C                 decimal place accuracy in the summation below.
                  SUM = C(NTERMS)
                  DO 60 K = NTERMS - 1, 1, -1
                     SUM = A*SUM + C(K)
   60             CONTINUE
                  GA = -SUM/(ONE+A*SUM)
                  Y = A*ALX
                  SUM = ONE
                  TERM = ONE
                  K = 1
   80             K = K + 1
                  IF (K.GT.MAXIT) GO TO 120
                  TERM = Y*TERM/K
                  SUM = SUM + TERM
                  IF (ABS(TERM).GT.EPS1*SUM) GO TO 80
                  U = GA - SUM*ALX
               END IF
               PP = A*X
               QQ = AP1
               RR = A + THREE
               TERM = ONE
               SUM = ONE
               K = 1
  100          K = K + 1
               IF (K.GT.MAXIT) GO TO 120
               PP = PP + X
               QQ = QQ + RR
               RR = RR + TWO
               TERM = -PP*TERM/QQ
               SUM = SUM + TERM
               IF (ABS(TERM).GT.EPS1*SUM) GO TO 100
               V = (X**AP1)*SUM/AP1
               Q = U + V
               Q = A*Q*EXP(-ALGP1)
               P = ONE - Q
            END IF
         END IF
      END IF
      GO TO 140
  120 IERR = 3
      NREC = 1
      WRITE (REC,FMT=99997) MAXIT
  140 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, A.le.0.0 : A =',1P,D13.5)
99998 FORMAT (1X,'** On entry, X.lt.0.0 : X =',1P,D13.5)
99997 FORMAT (1X,'** Algorithm fails to terminate in ',I4,' iterations')
      END
