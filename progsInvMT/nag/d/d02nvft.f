      SUBROUTINE D02NVF(NEQMAX,NY2DIM,MAXORD,METHOD,PETZLD,CONST,TCRIT,
     *                  HMIN,HMAX,H0,MAXSTP,MXHNIL,NORM,RWORK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     SETUP ROUTINE FOR INTEGRATOR EMPLOYING B.D.F. METHOD
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NVF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H0, HMAX, HMIN, TCRIT
      INTEGER           IFAIL, MAXORD, MAXSTP, MXHNIL, NEQMAX, NY2DIM
      LOGICAL           PETZLD
      CHARACTER*1       METHOD, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  CONST(6), RWORK(50+4*NEQMAX)
C     .. Local Scalars ..
      INTEGER           IDEV, IERR, NMETH
      LOGICAL           REPORT
      CHARACTER*1       METHD1
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NVZ, E04UDU, X04AAF, X04BAF
C     .. Executable Statements ..
      IERR = 0
      REPORT = IFAIL .LE. 0
      CALL X04AAF(0,IDEV)
      IF (NEQMAX.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99998)
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      METHD1 = METHOD
C     ENSURE UPPER CASE
      CALL E04UDU(METHD1)
      IF (METHD1.EQ.'F') THEN
         NMETH = 20
      ELSE IF (METHD1.EQ.'N' .OR. METHD1.EQ.'D') THEN
         NMETH = 2
      ELSE
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99999) METHD1
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      CALL D02NVZ(NEQMAX,NY2DIM,MAXORD,NMETH,PETZLD,CONST,TCRIT,HMIN,
     *            HMAX,H0,MAXSTP,MXHNIL,NORM,RWORK,IERR,REPORT)
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
99999 FORMAT (' ** D02NVF - METHOD(1:1)(=',A1,') IN NOT ''N'' OR ''F''',
     *  ' OR ''D'' **')
99998 FORMAT (' ** D02NVF - NEQMAX(=',I16,') .LT. 1 ** ')
      END
