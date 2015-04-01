      SUBROUTINE G01EPF(N,IP,D,PDL,PDU,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Calculates upper and lower bounds for the significance of a
C     Durbin-Watson statistic.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01EPF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D, PDL, PDU
      INTEGER           IFAIL, IP, N
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, PI
      INTEGER           I, IERROR, NIP, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01JDY, G01JDZ, X01AAF
      INTEGER           P01ABF
      EXTERNAL          G01JDY, G01JDZ, X01AAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
C
C
      IF (N.LE.IP) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N, IP
      ELSE IF (IP.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) IP
      ELSE IF (D.LT.0.0D0) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99997) D
      ELSE
C
         C = 0.0D0
         PI = X01AAF(0.0D0)
         NIP = N - IP
C
         DO 20 I = 1, N - 1
            WORK(I) = 2.0D0*(1.0D0-COS(PI*DBLE(I)/DBLE(N))) - D
   20    CONTINUE
         IF (N.LE.60) THEN
            PDL = G01JDZ(NIP,WORK,C)
         ELSE
            PDL = G01JDY(NIP,WORK,C)
         END IF
         DO 40 I = 1, N - 1
            WORK(I) = 2.0D0*(1.0D0-COS(PI*DBLE(I+IP-1)/DBLE(N))) - D
   40    CONTINUE
         IF (N.LE.60) THEN
            PDU = G01JDZ(NIP,WORK,C)
         ELSE
            PDU = G01JDY(NIP,WORK,C)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.le.IP : N = ',I16,' IP = ',I16)
99998 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99997 FORMAT (' ** On entry, D.lt.0.0 : D = ',D13.5)
      END
