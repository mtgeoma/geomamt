      SUBROUTINE D02NTF(NEQ,NEQMAX,JCEVAL,ML,MU,NWKJAC,NJCPVT,RWORK,
     *                  IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     SETUP ROUTINE FOR BANDED LINEAR ALGEBRA OPTION
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NTF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, ML, MU, NEQ, NEQMAX, NJCPVT, NWKJAC
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
      CHARACTER*80      REC(3)
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
      LAOPTN = 'BANDED'
      IERR = 0
      REPORT = IFAIL .LE. 0
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
            WRITE (REC(1),FMT=99989) JCEVL1
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
      IF (ML.LT.0) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99991) ML
            CALL X04BAF(IDEV,REC(1))
         END IF
      ELSE IF (ML.GT.NEQ-1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99992) ML, NEQ - 1
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (MU.LT.0) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(2),FMT=99993) MU
            CALL X04BAF(IDEV,REC(2))
         END IF
      ELSE IF (MU.GT.NEQ-1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(2),FMT=99990) MU, NEQ - 1
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
      IF (NJCPVT.LT.NEQMAX) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(2),FMT=99996) NJCPVT, NEQMAX
            CALL X04BAF(IDEV,REC(2))
         END IF
      END IF
      LENREQ = (2*ML+MU+1)*NEQMAX
      IF (NWKJAC.LT.LENREQ) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(3),FMT=99995) NWKJAC
            CALL X04BAF(IDEV,REC(3))
            WRITE (REC(3),FMT=99994) LENREQ
            CALL X04BAF(IDEV,REC(3))
         END IF
      END IF
C
      IF (IERR.EQ.0) THEN
         RWORK(34) = 2.0D0
         RWORK(35) = 1.0D0
         RWORK(36) = DBLE(ML)
         RWORK(37) = DBLE(MU)
         RWORK(46) = DBLE(NWKJAC)
         RWORK(47) = DBLE(NJCPVT)
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
99999 FORMAT (' ** D02NTF - NEQ(=',I16,') .LT. 1 **')
99998 FORMAT (' ** D02NTF - NEQMAX(=',I16,') .LT. 1 **')
99997 FORMAT (' ** D02NTF - NEQ(=',I16,') .GT. NEQMAX(=',I16,') **')
99996 FORMAT (' ** D02NTF - NJCPVT(=',I16,') .LT. NEQMAX(=',I16,') **')
99995 FORMAT (' ** D02NTF - NWKJAC(=',I16,') .LT. LENGTH REQUIRED')
99994 FORMAT ('             WHICH IS (2*ML+MU+1)*NEQMAX(=',I16,') **')
99993 FORMAT (' ** D02NTF - ML(=',I16,') .LT. 0 **')
99992 FORMAT (' ** D02NTF - ML(=',I16,') .GT. NEQ-1(=',I16,') **')
99991 FORMAT (' ** D02NTF - MU(=',I16,') .LT. 0 **')
99990 FORMAT (' ** D02NTF - MU(=',I16,') .GT. NEQ-1(=',I16,') **')
99989 FORMAT (' ** D02NTF - JCEVAL(1:1)(=',A1,') IS NOT ''A'' OR ''N''',
     *  ' OR ''D'' **')
      END
