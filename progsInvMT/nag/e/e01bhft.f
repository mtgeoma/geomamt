      SUBROUTINE E01BHF(N,X,F,D,A,B,PINT,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Evaluates the definite integral of a piecewise cubic Hermite
C     interpolant over the interval (a,b).
C
C     E01BHF is a driver for E01BHZ (derived from PCHIP routine PCHIA),
C     specialised for the case INCFD = 1.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01BHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, PINT
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), F(N), X(N)
C     .. Local Scalars ..
      INTEGER           I, IERR, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  E01BHZ
      INTEGER           P01ABF
      EXTERNAL          E01BHZ, P01ABF
C     .. Executable Statements ..
      IERR = 0
      IF (N.LT.2) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) N
         GO TO 40
      ELSE
         DO 20 I = 2, N
            IF (X(I).LE.X(I-1)) THEN
               IERR = 2
               NREC = 2
               WRITE (REC,FMT=99998) I, X(I-1), X(I)
               GO TO 40
            END IF
   20    CONTINUE
      END IF
C
      PINT = E01BHZ(N,X,F,D,1,A,B,IERR)
C
C     If A or B was out of range, signal that result is suspect.
      IF (IERR.LT.0) THEN
         IERR = 3
         NREC = 3
         WRITE (REC,FMT=99997) A, B
      END IF
C
   40 IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
99999 FORMAT (' ** On entry, N .lt. 2: N =',I16,'.')
99998 FORMAT (' ** On entry, X(R-1) .ge. X(R) for R =',I10,':',/'    X',
     *       '(R-1), X(R) =',1P,2D13.5,' .')
99997 FORMAT (' ** Warning - either A or B is outside the range X(1) .',
     *       '. X(N).',/4X,'The result has been computed by extrapolat',
     *       'ion and is unreliable.',/'    A =',1P,D13.5,'  B =',D13.5)
      END
