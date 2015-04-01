      SUBROUTINE E01BEF(N,X,F,D,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes a monotonicity-preserving piecewise cubic Hermite
C     interpolant to a set of data points.
C
C     E01BEF is a driver for E01BEZ (derived from PCHIP routine PCHIM),
C     specialiSed for the case INCFD = 1.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01BEF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), F(N), X(N)
C     .. Local Scalars ..
      INTEGER           I, IDUMM, IERR, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01BEZ
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
      CALL E01BEZ(N,X,F,D,1,IDUMM)
C
   40 IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
99999 FORMAT (1X,'** On entry, N .lt. 2: N =',I16,'.')
99998 FORMAT (1X,'** On entry, X(R-1) .ge. X(R) for R =',I10,':',/4X,
     *       'X(R-1), X(R) =',1P,2D13.5,' .')
      END
