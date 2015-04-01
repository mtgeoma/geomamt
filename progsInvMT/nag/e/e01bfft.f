      SUBROUTINE E01BFF(N,X,F,D,M,PX,PF,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Evaluates a piecewise cubic Hermite interpolant at a set of
C     points.
C
C     E01BFF is a driver for E01BFZ (derived from PCHIP routine PCHFE),
C     specialised for the case INCFD = 1.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01BFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), F(N), PF(M), PX(M), X(N)
C     .. Local Scalars ..
      INTEGER           I, IERR, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01BFZ
C     .. Executable Statements ..
C     Validity-check arguments.
      IERR = 0
      IF (N.LT.2) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) N
         GO TO 40
      ELSE IF (M.LT.1) THEN
         IERR = 3
         NREC = 1
         WRITE (REC,FMT=99997) M
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
      CALL E01BFZ(N,X,F,D,1,M,PX,PF,IERR)
C
C     If IERR .lt. 0, signal that extrapolation was performed at one or
C     more points.
C
      IF (IERR.LT.0) THEN
         IERR = 4
         NREC = 2
         WRITE (REC,FMT=99996)
      END IF
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
99997 FORMAT (1X,'** On entry, M .lt. 1: M =',I16,'.')
99996 FORMAT (1X,'** Warning - some points in array PX lie outside the',
     *       ' range X(1) .. X(N).',/4X,'Values at these points are un',
     *       'reliable because computed by extrapolation.')
      END
