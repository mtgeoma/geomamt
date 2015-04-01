      DOUBLE PRECISION FUNCTION G01GEF(X,A,B,RLAMDA,TOL,MAXIT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Based on algorithm AS226 APPL. STATIST. (1987) VOL. 36, NO. 2
C
C     Returns the cumulative probability of X for the non - central
C     Beta Distribution with parameters A, B and noncentrality RLAMDA.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, HALF, ONE, BIG, TEN
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,
     *                                 BIG=1.0D6,TEN=10.0D0)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01GEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, RLAMDA, TOL, X
      INTEGER                          IFAIL, MAXIT
C     .. Local Scalars ..
      DOUBLE PRECISION                 A1, AX, C, ERRBD, GX, P, PDF,
     *                                 PREC, Q, SUMQ, TGX, UFLO, XJ
      INTEGER                          IERR, IFAIL2, IFAULT
      LOGICAL                          IND1, IND2
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 S14ABF, X02AJF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         S14ABF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         G01BJW, G01EEF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, INT, LOG, MAX, MIN
C     .. Executable Statements ..
C
C      Define UFLO for underflow traps.
C
      UFLO = LOG(X02AMF())
C
      G01GEF = ZERO
      IERR = 1
      IF (A.LE.ZERO .OR. B.LE.ZERO) THEN
         WRITE (REC,FMT=99998) A, B
      ELSE IF (A.GT.BIG .OR. B.GT.BIG) THEN
         WRITE (REC,FMT=99993) A, B
      ELSE IF (RLAMDA.LT.ZERO .OR. RLAMDA.GT.-2.0D0*UFLO) THEN
         WRITE (REC,FMT=99997) RLAMDA
      ELSE IF (X.LT.ZERO .OR. X.GT.ONE) THEN
         WRITE (REC,FMT=99996) X
      ELSE IF (MAXIT.LT.1) THEN
         WRITE (REC,FMT=99992) MAXIT
      ELSE
         IERR = 0
      END IF
      IF (IERR.EQ.1) GO TO 40
      PREC = X02AJF()*10.0D0
      IF (TOL.GT.PREC .AND. TOL.LT.ONE) PREC = TOL
      G01GEF = X
C
C      Evaluate the central case
C
      IF (RLAMDA.EQ.ZERO) THEN
         IFAIL2 = 1
         CALL G01EEF(X,A,B,PREC,P,Q,PDF,IFAIL2)
         G01GEF = P
         IF (IFAIL2.EQ.4) THEN
            IERR = 3
            WRITE (REC,FMT=99995)
         ELSE IF (IFAIL2.EQ.3) THEN
            IERR = 4
            WRITE (REC,FMT=99994)
         END IF
      ELSE IF (X.NE.ZERO .AND. X.NE.ONE) THEN
         C = RLAMDA*HALF
C
C        Initialize the series...
C
         IND1 = (.FALSE.)
         IND2 = (.FALSE.)
         IF (A.GE.TEN .OR. B.GE.TEN) THEN
            CALL G01BJW(X,A,B,TGX)
            TGX = TGX + LOG(X) + LOG(ONE-X) - LOG(A)
         ELSE
            IFAULT = 0
            TGX = S14ABF(A+B,IFAULT) - S14ABF(A,IFAULT) - S14ABF(B,
     *            IFAULT) + A*LOG(X) + B*LOG(ONE-X) - LOG(A)
         END IF
         IFAIL2 = 1
         CALL G01EEF(X,A,B,PREC,P,Q,PDF,IFAIL2)
         IF (IFAIL2.EQ.3) IND2 = (.TRUE.)
C
C        Check for extreme cases
C
         IF (TGX.LT.(UFLO)) THEN
            IF (P.EQ.ZERO) THEN
               WRITE (REC,FMT=99995)
               IERR = 3
               G01GEF = ZERO
               GO TO 40
            ELSE IF (P.EQ.ONE) THEN
               A1 = A + MAXIT
               CALL G01EEF(X,A1,B,PREC,P,Q,PDF,IFAIL2)
               IF (P.EQ.ONE) THEN
                  WRITE (REC,FMT=99995)
                  IERR = 3
                  G01GEF = ONE
                  GO TO 40
               ELSE
C
C                 Initialize GX
C
                  GX = ZERO
                  IND1 = (.TRUE.)
               END IF
            ELSE
               GX = ZERO
               IND1 = (.TRUE.)
            END IF
         ELSE
            GX = EXP(TGX)
         END IF
         Q = EXP(-C)
         XJ = ZERO
         AX = Q*P
         SUMQ = ONE - Q
         G01GEF = AX
   20    XJ = XJ + ONE
         P = P - GX
         IF (IND1) THEN
            A1 = A + XJ
            IF (B.GT.ONE) THEN
               CALL G01BJW(X,A1,B,TGX)
               TGX = TGX + LOG(X) + LOG(ONE-X) - LOG(A1)
            ELSE
               CALL G01BJW(X,A1,B+ONE,TGX)
               TGX = TGX + LOG(X) + LOG(B/(A1*(A1+B)))
            END IF
            IF (TGX.LT.UFLO) THEN
               GX = ZERO
            ELSE
               GX = EXP(TGX)
               IND1 = (.FALSE.)
            END IF
         ELSE
            GX = X*(A+B+XJ-ONE)*GX/(A+XJ)
         END IF
         Q = Q*C/XJ
         SUMQ = SUMQ - Q
         AX = P*Q
         G01GEF = G01GEF + AX
C
C        Check for convergence and act accordingly...
C
         ERRBD = (P-GX)*SUMQ
         IF ((INT(XJ).LT.MAXIT) .AND. (ABS(ERRBD).GT.PREC)) GO TO 20
         IF (ERRBD.GT.PREC) THEN
            WRITE (REC,FMT=99999) MAXIT
            IERR = 2
         ELSE IF (IND2) THEN
            IERR = 4
            WRITE (REC,FMT=99994)
         END IF
         G01GEF = MAX(MIN(G01GEF,ONE),ZERO)
      END IF
   40 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** The solution has failed to converge in ',I16,' it',
     *       'erations.')
99998 FORMAT (1X,'** On entry A.le.0.0 or B.le.0.0: A = ',D13.5,' B = ',
     *       D13.5)
99997 FORMAT (1X,'** On entry, RLAMDA.lt.0.0 or RLAMDA too large: RLAM',
     *       'DA = ',1P,D13.5)
99996 FORMAT (1X,'** On entry, X.lt.0.0 OR .X.gt.1.0: X = ',1P,D13.5)
99995 FORMAT (1X,'** Probability too close to 0.0 or 1.0.')
99994 FORMAT (1X,'** Required accuracy not achieved in calculating cen',
     *       'tral beta probabilities')
99993 FORMAT (1X,'** On entry A or B too big: A = ',D13.5,' B = ',D13.5)
99992 FORMAT (1X,'** On entry, MAXIT.lt.1: MAXIT = ',I13)
      END
