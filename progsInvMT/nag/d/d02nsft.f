      SUBROUTINE D02NSF(NEQ,NEQMAX,JCEVAL,NWKJAC,RWORK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     SETUP ROUTINE FOR FULL LINEAR ALGEBRA OPTION
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NSF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NEQ, NEQMAX, NWKJAC
      CHARACTER*1       JCEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(50+4*NEQMAX)
C     .. Scalars in Common ..
      CHARACTER*6       LAOPTN
C     .. Local Scalars ..
      INTEGER           IDEV, IERR, LENREQ
      LOGICAL           REPORT
      CHARACTER*1       JCEVL1
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04UDU, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Common blocks ..
      COMMON            /AD02NC/LAOPTN
C     .. Save statement ..
      SAVE              /AD02NC/
C     .. Executable Statements ..
      LAOPTN = 'FULLOP'
C     REPORT IS USED TO SIMULATE QUIET/NOISY FAILURES IN P01ABF
      REPORT = IFAIL .LE. 0
      IERR = 0
      CALL X04AAF(0,IDEV)
      JCEVL1 = JCEVAL
C     ENSURE UPPER CASE
      CALL E04UDU(JCEVL1)
      IF (JCEVL1.EQ.'A') THEN
         RWORK(38) = 1.0D0
      ELSE IF (JCEVL1.EQ.'N' .OR. JCEVL1.EQ.'D') THEN
         RWORK(38) = 0.0D0
      ELSE
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99994) JCEVL1
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
C
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
            WRITE (REC(2),FMT=99998) NEQMAX
            CALL X04BAF(IDEV,REC(2))
         END IF
      END IF
C
      IF (NEQ.GT.NEQMAX) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99997) NEQ, NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      LENREQ = NEQMAX*(NEQMAX+1)
      IF (NWKJAC.LT.LENREQ) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(2),FMT=99996) NWKJAC
            CALL X04BAF(IDEV,REC(2))
            WRITE (REC(2),FMT=99995) LENREQ
            CALL X04BAF(IDEV,REC(2))
         END IF
      END IF
C
      IF (IERR.EQ.0) THEN
         RWORK(34) = 1.0D0
         RWORK(35) = 1.0D0
         RWORK(46) = DBLE(NWKJAC)
         RWORK(47) = 1.0D0
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
99999 FORMAT (' ** D02NSF - NEQ(=',I16,') .LT. 1 **')
99998 FORMAT (' ** D02NSF - NEQMAX(=',I16,') .LT. 1 **')
99997 FORMAT (' ** D02NSF - NEQ(=',I16,') .GT. NEQMAX(=',I16,') **')
99996 FORMAT (' ** D02NSF - NWKJAC(=',I16,') .LT. REQUIRED LENGTH ')
99995 FORMAT ('             WHICH IS NEQMAX*(NEQMAX+1) (=',I16,') **')
99994 FORMAT (' ** D02NSF - JCEVAL(1:1)(=',A1,') IS NOT ''N'' OR ''A''',
     *  ' OR ''D'' **')
      END
