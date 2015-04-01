      DOUBLE PRECISION FUNCTION G01EDF(TAIL,F,DF1,DF2,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Computes the lower/upper tail probability of F for the
C     F-distribution with DF1 and DF2 degrees of freedom.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, TWO, ONE, BIG
      PARAMETER                        (ZERO=0.0D0,TWO=2.0D0,ONE=1.0D0,
     *                                 BIG=1.0D5)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01EDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF1, DF2, F
      INTEGER                          IFAIL
      CHARACTER                        TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 P, PDF, Q, TOL, X
      INTEGER                          IERR, IFAIL2
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01ECF, S15ABF
      INTEGER                          P01ABF
      EXTERNAL                         G01ECF, S15ABF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         G01EEF
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Executable Statements ..
C
      G01EDF = ZERO
C
C     Check the input parameters.
C
      IERR = 0
      IF (DF1.LE.ZERO .OR. DF2.LE.ZERO) THEN
         IERR = 3
         WRITE (REC,FMT=99999) DF1, DF2
      ELSE IF (F.LT.ZERO) THEN
         IERR = 2
         WRITE (REC,FMT=99998) F
      ELSE IF (TAIL.NE.'L' .AND. TAIL.NE.'U' .AND. TAIL.NE.'l' .AND.
     *         TAIL.NE.'u') THEN
         IERR = 1
         WRITE (REC,FMT=99997) TAIL
      ELSE IF (DF1.GT.BIG .AND. DF2.GT.BIG) THEN
C
C        Use normal approximation
C
         P = 2.0D0/(9.0D0*DF1)
         Q = 2.0D0/(9.0D0*DF2)
         X = F**(1.0D0/3.0D0)
         IFAIL2 = 1
         IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
            G01EDF = S15ABF(-(X*(ONE-Q)-ONE+P)/SQRT(P+X*X*Q),IFAIL2)
         ELSE
            G01EDF = S15ABF((X*(ONE-Q)-ONE+P)/SQRT(P+X*X*Q),IFAIL2)
         END IF
C
C        Use Chisquare approximation
C
      ELSE IF (DF2.GT.BIG) THEN
         IFAIL2 = 1
         G01EDF = G01ECF(TAIL,DF1*F,DF1,IFAIL2)
      ELSE IF (DF1.GT.BIG) THEN
         IFAIL2 = 1
         IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
            G01EDF = G01ECF('L',DF2/F,DF2,IFAIL2)
         ELSE
            G01EDF = G01ECF('U',DF2/F,DF2,IFAIL2)
         END IF
      ELSE
         TOL = 0.5D-5
         X = DF1*F/(DF1*F+DF2)
         IFAIL2 = 1
C
C        Use a transformation of a BETA.
C
         CALL G01EEF(X,DF1/TWO,DF2/TWO,TOL,P,Q,PDF,IFAIL2)
         IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
            G01EDF = P
         ELSE
            G01EDF = Q
         END IF
         IF (IFAIL2.EQ.4) THEN
            IERR = 4
            WRITE (REC,FMT=99996)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, DF1.le.0.0 or DF2.le.0.0: DF1 = ',D13.5,
     *       ' DF2 = ',D13.5)
99998 FORMAT (1X,'** On entry, F.lt.0.0: F = ',D13.5)
99997 FORMAT (1X,'** On entry TAIL is not valid: TAIL = ',A1)
99996 FORMAT (1X,'** F is too far out into the tails for accurate calc',
     *       'ulation')
      END
