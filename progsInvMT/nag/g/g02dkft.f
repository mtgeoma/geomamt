      SUBROUTINE G02DKF(IP,ICONST,P,C,LDC,B,RSS,IDF,SE,COV,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES LEAST SQUARES ESTIMATES FOR GIVEN SET OF LINEAR
C     CONSTRAINTS GIVEN THE SVD SOLUTION FROM G02DAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02DKF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS
      INTEGER           ICONST, IDF, IFAIL, IP, LDC
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), C(LDC,ICONST), COV((IP*IP+IP)/2),
     *                  P(IP*(IP+2)), SE(IP), WK(2*IP*IP+IP*ICONST+2*
     *                  ICONST*ICONST+4*ICONST)
C     .. Local Scalars ..
      DOUBLE PRECISION  S
      INTEGER           IERROR, IND, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02DKZ
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (IP.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) IP
      ELSE IF (ICONST.GE.IP) THEN
         WRITE (P01REC(1),FMT=99998) ICONST, IP
      ELSE IF (ICONST.LE.0) THEN
         WRITE (P01REC(1),FMT=99995) ICONST
      ELSE IF (LDC.LT.IP) THEN
         WRITE (P01REC(1),FMT=99996) LDC, IP
      ELSE IF (RSS.LE.0.0D0) THEN
         WRITE (P01REC(1),FMT=99997) RSS
      ELSE IF (IDF.LE.0) THEN
         WRITE (P01REC(1),FMT=99993) IDF
      ELSE
         IERROR = 0
         S = RSS/IDF
         CALL G02DKZ(IP,ICONST,P,0,C,LDC,B,S,SE,COV,WK,IND)
         IF (IND.NE.0) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99994)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99998 FORMAT (' ** On entry, ICONST.ge.IP : ICONST = ',I16,' IP = ',I16)
99997 FORMAT (' ** On entry, RSS.le.0.0 : RSS = ',D13.5)
99996 FORMAT (' ** On entry, LDC.lt.IP : LDC = ',I16,' IP = ',I16)
99995 FORMAT (' ** On entry, ICONST.le.0 : ICONST =',I16)
99994 FORMAT (' ** C does not give a model of full rank')
99993 FORMAT (' ** On entry, IDF.le.0 : IDF =',I16)
      END
