      SUBROUTINE D02NYF(NEQ,NEQMAX,HU,H,TCUR,TOLSF,RWORK,NST,NRE,NJE,
     *                  NQU,NQ,NITER,IMXER,ALGEQU,INFORM,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-531 (FEB 1987).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NYF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HU, TCUR, TOLSF
      INTEGER           IFAIL, IMXER, NEQ, NEQMAX, NITER, NJE, NQ, NQU,
     *                  NRE, NST
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(50+4*NEQMAX)
      INTEGER           INFORM(23)
      LOGICAL           ALGEQU(NEQ)
C     .. Local Scalars ..
      INTEGER           I, IALGQB, IDEV, IERR
      LOGICAL           REPORT
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Executable Statements ..
      IERR = 0
      REPORT = IFAIL .LE. 0
      CALL X04AAF(0,IDEV)
      IF (NEQ.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99999) NEQ
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (NEQMAX.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99998) NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (IERR.EQ.0 .AND. NEQ.GT.NEQMAX) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99997) NEQ, NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         HU = RWORK(11)
         H = RWORK(12)
         TCUR = RWORK(13)
         TOLSF = RWORK(14)
C
         NST = INFORM(1)
         NRE = INFORM(2)
         NJE = INFORM(3)
         NQU = INFORM(4)
         NQ = INFORM(5)
         NITER = INFORM(6)
         IMXER = INFORM(8)
         IALGQB = 50 + 3*NEQMAX
         DO 20 I = 1, NEQ
            ALGEQU(I) = RWORK(IALGQB+I) .EQ. 0.0D0
   20    CONTINUE
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
99999 FORMAT (' ** D02NYF - NEQ(=',I16,') .LT. 1 **')
99998 FORMAT (' ** D02NYF - NEQMAX(=',I16,') .LT. 1 **')
99997 FORMAT (' ** D02NYF - NEQ(=',I16,') .GT. NEQMAX(=',I16,') **')
      END
