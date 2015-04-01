      DOUBLE PRECISION FUNCTION G01FEF(P,A,B,TOL,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-918 (APR 1991).
C
C     G01FEF RETURNS THE DEVIATE ASSOCIATED
C     WITH THE LOWER TAIL PROBABILITY P FROM
C     THE BETA DISTRIBUTION OF THE FIRST KIND
C     WITH PARAMETERS A AND B.
C
C     BASED ON ALGORITHM AS 109 APPL. STATIST. (1977), VOL.26, NO.1
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FEF')
      DOUBLE PRECISION                 ZERO, ONE, TEN, BIG
      PARAMETER                        (ZERO=0.0D0,ONE=1.0D0,TEN=10.0D0,
     *                                 BIG=1.0D6)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, P, TOL
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 ADJ, BINV, BLOG, EPS, G, H,
     *                                 POWER, PR, R, S, SEPS, SQ, T, TX,
     *                                 UFLO, W, Y, YP, YY, YYY
      INTEGER                          IERR, IFA, ITER, MAXIT
      LOGICAL                          NACC
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 S14ABF, X02AJF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         S14ABF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         G01EEF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, LOG, MAX, MIN, SIGN,
     *                                 SQRT
C     .. Executable Statements ..
C
C     CHECK PARAMETER RANGES
C
      G01FEF = ZERO
      MAXIT = 100
      IERR = 0
C
      IF (P.LT.ZERO .OR. P.GT.ONE) THEN
         IERR = 1
         WRITE (REC,FMT=99999) P
         GO TO 200
      ELSE IF (A.LE.ZERO .OR. B.LE.ZERO) THEN
         IERR = 2
         WRITE (REC,FMT=99998) A, B
         GO TO 200
      ELSE IF (A.GT.BIG .OR. B.GT.BIG) THEN
         IERR = 2
         WRITE (REC,FMT=99996) A, B
         GO TO 200
      END IF
      NACC = .FALSE.
      G01FEF = P
      EPS = X02AJF()
      SEPS = SQRT(EPS)
      EPS = EPS*TEN
      IF (TOL.GT.EPS .AND. TOL.LT.ONE) EPS = TOL
      UFLO = LOG(X02AMF())
      IF (P.NE.ZERO .AND. P.NE.ONE) THEN
C
C        CALCULATE LN(BETA(A,B))
C
         IFA = 1
         BLOG = (S14ABF(A,IFA)+S14ABF(B,IFA)) - S14ABF(A+B,IFA)
C
C        FIND INITIAL APPROXIMATION
C
         R = SQRT(-LOG(P*P))
         Y = R - (2.30753D0+0.27061D0*R)/(ONE+(0.99229D0+0.04481D0*R)*R)
         IF (A.GT.ONE .AND. B.GT.ONE) GO TO 60
         R = B + B
         T = ONE/(9.0D0*B)
         T = R*(ONE-T+Y*SQRT(T))**3
         IF (T.LE.ZERO) GO TO 20
         T = (4.0D0*A+R-2.0D0)/T
         IF (T.LE.ONE) GO TO 40
         BINV = ONE - 2.0D0/(T+ONE)
         GO TO 80
   20    POWER = (LOG((ONE-P)*B)+BLOG)/B
         IF (POWER.GT.ZERO) THEN
            BINV = SEPS
         ELSE IF (POWER.LT.UFLO) THEN
            BINV = ONE - SEPS
         ELSE
            BINV = ONE - EXP(POWER)
         END IF
         GO TO 80
   40    POWER = (LOG(P*A)+BLOG)/A
         IF (POWER.LT.UFLO) THEN
            BINV = SEPS
         ELSE
            BINV = EXP(POWER)
         END IF
         GO TO 80
   60    R = (Y*Y-3.0D0)/6.0D0
         S = ONE/(A+A-ONE)
         T = ONE/(B+B-ONE)
         H = 2.0D0/(S+T)
         W = Y*SQRT(H+R)/H - (T-S)*(R+5.0D0/6.0D0-2.0D0/(3.0D0*H))
         BINV = A/(A+B*EXP(W+W))
C
C        CONVERGE USING MODIFIED NEWTON-RAPHSON METHOD
C
   80    YP = ZERO
         SQ = ONE
         PR = ONE
         BINV = MIN(MAX(SEPS,BINV),ONE-SEPS)
         IFA = 1
         DO 160 ITER = 1, MAXIT
            IFA = 1
            CALL G01EEF(BINV,A,B,EPS,Y,YY,YYY,IFA)
            IF (IFA.NE.0) NACC = .TRUE.
            Y = (Y-P)*EXP(BLOG-(A-ONE)*LOG(BINV)-(B-ONE)*LOG(ONE-BINV))
            IF (SIGN(ONE,Y).NE.SIGN(ONE,YP)) PR = SQ
            G = ONE
  100       ADJ = G*Y
            SQ = ABS(ADJ)
            IF (SQ.GE.PR) GO TO 120
            TX = BINV - ADJ
            IF (TX.GE.ZERO .AND. TX.LE.ONE) GO TO 140
  120       G = G/3.0D0
            GO TO 100
  140       SEPS = EPS*BINV
            IF (PR.LE.SEPS .OR. ABS(Y).LE.SEPS) GO TO 180
            IF (TX.EQ.ZERO .OR. TX.EQ.ONE) GO TO 120
            IF (TX.EQ.BINV) GO TO 180
            BINV = TX
            YP = Y
  160    CONTINUE
         IERR = 3
         WRITE (REC,FMT=99997)
  180    G01FEF = BINV
         IF (IERR.EQ.0 .AND. NACC) THEN
            IERR = 4
            WRITE (REC,FMT=99995)
         END IF
      END IF
  200 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry P.lt.0.0 .or. P.gt.1.0: P = ',D13.5)
99998 FORMAT (1X,'** On entry either A or B .le. 0.0: A = ',D13.5,' B ',
     *       '= ',D13.5)
99997 FORMAT (1X,'** The solution failed to converge.')
99996 FORMAT (1X,'** On entry A and/or B is too large: A = ',D13.5,' B',
     *       ' = ',D13.5)
99995 FORMAT (1X,'** Requested accuracy not achieved')
      END
