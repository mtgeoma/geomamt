      DOUBLE PRECISION FUNCTION G01FFF(P,A,B,TOL,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-933 (APR 1991).
C
C     Based on:
C     Algorithm as 91 Appl. Statist. (1975) Vol.24, P.385.
C     Computes the deviate associated with the lower tail
C     probability P, of the Gamma distribution with shape
C     parameter A and scale parameter B.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, HALF, ONE, TWO, THREE, SIX,
     *                                 FIFTY, BIG
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,
     *                                 TWO=2.0D0,THREE=3.0D0,SIX=6.0D0,
     *                                 FIFTY=50.0D0,BIG=1.0D6)
      DOUBLE PRECISION                 C1, C2, C3, C4, C5, C6, C7, C8,
     *                                 C9, C10, C11, C12, C13, C14, C15,
     *                                 C16, C17, C18, C19, C20, C21,
     *                                 C22, C23, C24, C25, C26, C27,
     *                                 C28, C29, C30, C31, C32, C33,
     *                                 C34, C35, C36, C37, C38
      PARAMETER                        (C1=0.01D0,C2=0.111111D0,
     *                                 C3=0.16D0,C4=0.4D0,C5=0.62D0,
     *                                 C6=4.4D0,C7=4.67D0,C8=6.66D0,
     *                                 C9=6.73D0,C10=13.32D0,C11=60.0D0,
     *                                 C12=70.0D0,C13=84.0D0,
     *                                 C14=105.0D0,C15=120.0D0,
     *                                 C16=127.0D0,C17=140.0D0,
     *                                 C18=175.0D0,C19=210.0D0,
     *                                 C20=252.0D0,C21=264.0D0,
     *                                 C22=294.0D0,C23=346.0D0,
     *                                 C24=420.0D0,C25=462.0D0,
     *                                 C26=606.0D0,C27=672.0D0,
     *                                 C28=707.0D0,C29=735.0D0,
     *                                 C30=889.0D0,C31=932.0D0,
     *                                 C32=966.0D0,C33=1141.0D0,
     *                                 C34=1182.0D0,C35=1278.0D0,
     *                                 C36=1740.0D0,C37=2520.0D0,
     *                                 C38=5040.0D0)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, P, TOL
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 BB, C, D, DD, EPS, G, GAM,
     *                                 GAMLOG, P1, P2, Q, Q2, S1, S2,
     *                                 S3, S4, S5, S6, T, UFLOW, X
      INTEGER                          IERROR, IFAIL2, IFAULT, ITER,
     *                                 MAXIT
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01CEF, S14ABF, X02AJF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         G01CEF, S14ABF, X02AJF, X02AMF,
     *                                 P01ABF
C     .. External Subroutines ..
      EXTERNAL                         S14BAF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, LOG, MAX, SQRT
C     .. Executable Statements ..
C
C     Test arguments and initialize.
C
      G01FFF = ZERO
      IF (P.LT.ZERO .OR. P.GE.ONE) THEN
         WRITE (REC,FMT=99999) P
         IERROR = 1
      ELSE IF (A.LE.ZERO) THEN
         WRITE (REC,FMT=99998) A
         IERROR = 2
      ELSE IF (A.GT.BIG) THEN
         WRITE (REC,FMT=99993) A
         IERROR = 2
      ELSE IF (B.LE.ZERO) THEN
         WRITE (REC,FMT=99994) B
         IERROR = 2
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0 .AND. P.NE.0.0D0) THEN
         EPS = MAX(X02AJF(),1.0D-18)*FIFTY
         IF (TOL.GT.EPS .AND. TOL.LT.ONE) EPS = TOL
         MAXIT = 100
         UFLOW = LOG(X02AMF())
         DD = LOG(TWO)
         IFAULT = 1
         G = S14ABF(A,IFAULT)
         C = A - ONE
         IF (A.LT.-C5*LOG(P)) THEN
            GAMLOG = ONE/A*(LOG(P)+LOG(A)+G+A*DD)
            IF (GAMLOG.LT.UFLOW) THEN
               GAM = ZERO
               IERROR = 3
               WRITE (REC,FMT=99996)
               GO TO 100
            ELSE
               GAM = EXP(GAMLOG)
            END IF
            IF (GAM.LT.EPS) GO TO 80
            GO TO 40
         END IF
C
C        Find starting values
C
         IF (A.LE.C3) THEN
            GAM = C4
            D = LOG(ONE-P)
   20       Q = GAM
            P1 = ONE + GAM*(C7+GAM)
            P2 = GAM*(C9+GAM*(C8+GAM))
            T = -HALF + (C7+TWO*GAM)/P1 - (C9+GAM*(C10+THREE*GAM))/P2
            GAM = GAM - (ONE-EXP(D+G+HALF*GAM+C*DD)*P2/P1)/T
            IF (ABS(Q/GAM-ONE).GT.C1) GO TO 20
            GO TO 40
         END IF
C
C        Call G01CEF (P HAS BEEN CHECKED)
C
         IFAULT = 1
         X = G01CEF(P,IFAULT)
C
C        Starting approximation using Wilson and Hilferty estimate.
C
         P1 = C2/A
         GAM = TWO*A*(X*SQRT(P1)+ONE-P1)**3
C
C        Starting approximation for P tending to 1.0.
C
         IF (GAM.GT.C6*A+SIX) GAM = -TWO*(LOG(ONE-P)-C*LOG(HALF*GAM)+G)
C
C        Call S14BAF and calculate seven term taylor series.
C
   40    ITER = 0
   60    CONTINUE
         ITER = ITER + 1
         Q = GAM
         P1 = HALF*GAM
         IFAIL2 = 1
         CALL S14BAF(A,P1,EPS,P2,Q2,IFAIL2)
         IF (IFAIL2.NE.0) THEN
            IERROR = 5
            WRITE (REC,FMT=99995)
            GO TO 100
         END IF
         P2 = P - P2
         T = P2*EXP(A*DD+G+P1-C*LOG(GAM))
         BB = T/GAM
         D = HALF*T - BB*C
         S1 = (C19+D*(C17+D*(C14+D*(C13+D*(C12+C11*D)))))/C24
         S2 = (C24+D*(C29+D*(C32+D*(C33+C35*D))))/C37
         S3 = (C19+D*(C25+D*(C28+C31*D)))/C37
         S4 = (C20+D*(C27+C34*D)+C*(C22+D*(C30+C36*D)))/C38
         S5 = (C13+C21*D+C*(C18+C26*D))/C37
         S6 = (C15+C*(C23+C16*C))/C38
         GAM = GAM + T*(ONE+HALF*T*S1-BB*C*
     *         (S1-BB*(S2-BB*(S3-BB*(S4-BB*(S5-BB*S6))))))
         IF (ABS(Q/GAM-ONE).GT.EPS) THEN
            IF (ITER.GE.MAXIT) THEN
               IERROR = 4
               WRITE (REC,FMT=99997)
            ELSE
               GO TO 60
            END IF
         END IF
   80    G01FFF = HALF*GAM*B
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, P.lt.0.0 .OR. P.ge.1.0: P = ',D13.5)
99998 FORMAT (1X,'** On entry, A.le.0.0 : A = ',D13.5)
99997 FORMAT (1X,'** Solution fails to converge.')
99996 FORMAT (1X,'** P is too close to 0.0 or 1.0')
99995 FORMAT (1X,'** Convergence failure in calculating gamma integral.'
     *       )
99994 FORMAT (1X,'** On entry, B .le. 0.0 : B = ',D13.5)
99993 FORMAT (1X,'** On entry, A is too large: A = ',D13.5)
      END
