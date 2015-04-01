      SUBROUTINE G01EEF(X,A,B,TOL,P,Q,PDF,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-1031 (JUN 1993).
C
C        COMPUTES INCOMPLETE BETA FUNCTIONS P(a,b,x) and Q(a,b,x),
C        and also the PROBABILITY DENSITY FUNCTION.
C
C        BASED ON :-
C        ALGORITHM AS 63  APPL. STATIST. (1973) VOL.22, P.409
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01EEF')
      DOUBLE PRECISION  BIG, ONE, ZERO, TEN
      PARAMETER         (BIG=1.0D6,ONE=1.0D0,ZERO=0.0D0,TEN=10.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, P, PDF, Q, TOL, X
      INTEGER           IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, ACC, AI, APB, ATOL, BB, BETAP, CX, EPOWER,
     *                  EPS, POWER, RX, TEMP, TERM, UFLO, XX
      INTEGER           IERR, IFAULT, NS
      LOGICAL           INDEX
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S14ABF, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          S14ABF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01BJW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, MIN
C     .. Executable Statements ..
C
C        Test for admissibility of arguments
C
      IERR = 0
      IF (X.LT.ZERO .OR. X.GT.ONE) THEN
         IERR = 1
         WRITE (REC,FMT=99999) X
      ELSE IF (A.LE.ZERO .OR. B.LE.ZERO) THEN
         IERR = 2
         WRITE (REC,FMT=99998) A, B
      ELSE IF (A.GT.BIG .OR. B.GT.BIG) THEN
         IERR = 2
         WRITE (REC,FMT=99997) A, B
      ELSE
         EPS = X02AJF()
         ATOL = TOL
         IF (TOL.LT.10.0D0*EPS .OR. TOL.GE.ONE) ATOL = 10.0D0*EPS
         UFLO = LOG(X02AMF())
         IF (X.NE.ZERO .AND. X.NE.ONE) THEN
C
C           Change tail if necessary and determine S
C
            APB = A + B
            CX = ONE - X
            IF (A.GE.APB*X) THEN
               XX = X
               AA = A
               BB = B
               INDEX = .FALSE.
            ELSE
               XX = CX
               CX = X
               AA = B
               BB = A
               INDEX = .TRUE.
            END IF
            TERM = ONE
            AI = ONE
            BETAP = ONE
            NS = BB + CX*APB
C
C           Use Reduction Formulae of SOPER.
C
            RX = XX/CX
   20       CONTINUE
            TEMP = BB - AI
            IF (NS.EQ.0) RX = XX
   40       CONTINUE
            TEMP = TEMP*RX/(AA+AI)
            TERM = TERM*TEMP
            BETAP = BETAP + TERM
            IF (TERM.GT.ZERO) THEN
               ACC = (ONE-TEMP)*ATOL
               IF (BETAP*ACC.GE.TERM) THEN
                  GO TO 80
               ELSE IF (ACC.LT.EPS .AND. BETAP*EPS.GE.TERM) THEN
                  GO TO 60
               END IF
            ELSE IF (TERM.LT.ZERO .AND. ABS(TEMP).LT.ONE) THEN
               ACC = (ONE-TEMP+TERM)*BETAP*ATOL
               IF (ACC.GE.-TERM) THEN
                  GO TO 80
               ELSE IF (ABS(ACC).LT.EPS .AND. BETAP*EPS.GE.-TERM) THEN
                  GO TO 60
               END IF
            ELSE
               GO TO 80
            END IF
            AI = AI + ONE
            NS = NS - 1
            IF (NS.GE.0) THEN
               GO TO 20
            ELSE
               TEMP = APB
               APB = APB + ONE
               GO TO 40
            END IF
   60       IERR = 3
            WRITE (REC,FMT=99996)
C
C           Calculate result
C
   80       CONTINUE
            IF (AA.GT.TEN .AND. BB.GT.TEN) THEN
               CALL G01BJW(XX,AA,BB,POWER)
               POWER = POWER + LOG(XX)
            ELSE
               IFAULT = 0
               POWER = -S14ABF(AA,IFAULT) - S14ABF(BB,IFAULT) +
     *                 S14ABF(A+B,IFAULT) + AA*LOG(XX) + (BB-ONE)
     *                 *LOG(CX)
            END IF
            IF (POWER.LT.UFLO) THEN
               IERR = 4
               WRITE (REC,FMT=99995)
               IF ( .NOT. INDEX) THEN
                  P = ZERO
               ELSE
                  P = ONE
               END IF
               Q = ONE - P
               PDF = ZERO
               GO TO 100
            ELSE
               EPOWER = EXP(POWER)
               P = BETAP*EPOWER/AA
               P = MIN(MAX(P,0.0D0),1.0D0)
               PDF = EPOWER/XX
            END IF
            IF (INDEX) THEN
               Q = P
               P = ONE - Q
            ELSE
               Q = ONE - P
            END IF
         ELSE IF (X.EQ.ZERO) THEN
            P = ZERO
            Q = ONE
            PDF = ZERO
         ELSE
            P = ONE
            Q = ZERO
            PDF = ZERO
         END IF
      END IF
  100 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
C
99999 FORMAT (1X,'** On entry X.lt.0.0 or X.gt.1.0: X = ',D13.5)
99998 FORMAT (1X,'** On entry A.le.0.0 or B.le.0.0: A = ',D13.5,' B = ',
     *       D13.5)
99997 FORMAT (1X,'** On entry A or B is too large: A = ',D13.5,' B = ',
     *       D13.5)
99996 FORMAT (1X,'** Required accuracy not achieved')
99995 FORMAT (1X,'** Probability too close to 0.0 or 1.0.')
      END
