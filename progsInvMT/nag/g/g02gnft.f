      SUBROUTINE G02GNF(IP,IRANK,B,COV,V,LDV,F,EST,STAT,SESTAT,Z,TOL,WK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     GIVES MAXIMUM LIKELIHOOD ESTIMATE AND SE OF ESTIMABLE FUNCTION
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02GNF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SESTAT, STAT, TOL, Z
      INTEGER           IFAIL, IP, IRANK, LDV
      LOGICAL           EST
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV((IP*IP+IP)/2), F(IP), V(LDV,IP+7),
     *                  WK(IP)
C     .. Local Scalars ..
      INTEGER           IERROR, IND, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02DNZ
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (IP.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) IP
      ELSE IF (IRANK.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) IRANK
      ELSE IF (IP.LT.IRANK) THEN
         WRITE (P01REC(1),FMT=99997) IRANK, IP
      ELSE IF (LDV.LT.IP) THEN
         WRITE (P01REC(1),FMT=99994) LDV, IP
      ELSE
         IERROR = 0
         IF (IRANK.EQ.IP) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99996)
         END IF
         CALL G02DNZ(IP,IRANK,B,COV,V(1,8),LDV,F,EST,STAT,SESTAT,Z,TOL,
     *               WK,IND)
         IF (IND.EQ.1) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99995)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99998 FORMAT (' ** On entry, IRANK.lt.1 : IRANK = ',I16)
99997 FORMAT (' ** On entry, IRANK.gt.IP : IRANK = ',I16,' IP = ',I16)
99996 FORMAT (' ** On entry, IRANK=IP, ie model of full rank')
99995 FORMAT (' ** Standard error of statistic is 0.0')
99994 FORMAT (' ** On entry, LDV.lt.IP : LDV = ',I16,' IP = ',I16)
      END
