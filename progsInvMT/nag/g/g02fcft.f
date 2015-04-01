      SUBROUTINE G02FCF(N,IP,RES,D,PDL,PDU,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02FCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D, PDL, PDU
      INTEGER           IFAIL, IP, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RES(N), WORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  PREC, RN, RSS, SUMR, TEMP
      INTEGER           I, IERROR, IFAULT, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DDOT, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01EPF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
C
      IF (N.LE.IP) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N, IP
      ELSE IF (IP.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) IP
      ELSE
         SUMR = 0.0D0
         RN = DBLE(N)
         PREC = SQRT(X02AJF())
         DO 20 I = 1, N
            SUMR = SUMR + RES(I)
   20    CONTINUE
         TEMP = SUMR/RN
         IF (ABS(TEMP).GT.PREC) THEN
            IERROR = 2
            WRITE (P01REC,FMT=99997) TEMP
         ELSE
            RSS = DDOT(N,RES(1),1,RES(1),1)
            SUMR = 0.0D0
            DO 40 I = 1, N - 1
               SUMR = SUMR + (RES(I)-RES(I+1))**2
   40       CONTINUE
            IF (SUMR.LE.0) THEN
               IERROR = 3
               WRITE (P01REC(1),FMT=99996)
            ELSE
               D = SUMR/RSS
               IFAULT = 0
               CALL G01EPF(N,IP,D,PDL,PDU,WORK,IFAULT)
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.le.IP : N = ',I16,' IP = ',I16)
99998 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99997 FORMAT (' ** On entry, The mean of RES is not approximately 0.0,',
     *       ' mean = ',D13.5)
99996 FORMAT (' ** On entry all residuals are identical')
      END
