      SUBROUTINE D02NZF(NEQMAX,TCRIT,H,HMIN,HMAX,MAXSTP,MAXHNL,RWORK,
     *                  IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-532 (FEB 1987).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NZF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HMAX, HMIN, TCRIT
      INTEGER           IFAIL, MAXHNL, MAXSTP, NEQMAX
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(50+4*NEQMAX)
C     .. Local Scalars ..
      INTEGER           IDEV, IERR
      LOGICAL           REPORT
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      IERR = 0
      REPORT = IFAIL .LE. 0
      IF (NEQMAX.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            CALL X04AAF(0,IDEV)
            WRITE (REC(1),FMT=99999) NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         RWORK(1) = TCRIT
         RWORK(2) = DBLE(MAXSTP)
         IF (MAXHNL.GT.0) THEN
            RWORK(3) = DBLE(MAXHNL)
         ELSE
            RWORK(3) = 10.0D0
         END IF
         IF (HMAX.LT.0.0D0) THEN
            RWORK(6) = RWORK(18)
         ELSE
            RWORK(6) = HMAX
         END IF
         IF (HMIN.LT.0.0D0) THEN
            RWORK(7) = RWORK(17)
         ELSE
            RWORK(7) = HMIN
         END IF
         RWORK(8) = 2.0D0
         RWORK(9) = ABS(H)
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
99999 FORMAT (' ** D02NZF - NEQMAX(=',I16,') .LT. 1 **')
      END
