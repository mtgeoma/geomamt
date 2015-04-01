      SUBROUTINE G02ECF(MEAN,N,SIGSQ,TSS,NMOD,NTERMS,RSS,RSQ,CP,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES CP STATISTIC AND R-SQUARED VALUES
C     FROM RESIDUAL SUM OF SQUARES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02ECF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGSQ, TSS
      INTEGER           IFAIL, N, NMOD
      CHARACTER*1       MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  CP(NMOD), RSQ(NMOD), RSS(NMOD)
      INTEGER           NTERMS(NMOD)
C     .. Local Scalars ..
      DOUBLE PRECISION  CPI, RSSI
      INTEGER           I, IC, IERROR, INTP, IP, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (NMOD.LE.0) THEN
         WRITE (P01REC,FMT=99999) NMOD
      ELSE IF (SIGSQ.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99998) SIGSQ
      ELSE IF (TSS.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99997) TSS
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
            INTP = 0
         ELSE IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            INTP = 1
         ELSE
            IERROR = 1
            WRITE (P01REC,FMT=99996) MEAN
            GO TO 100
C
         END IF
C
C       CALCULATE CP AND R-SQUARED
C
         DO 20 I = 1, NMOD
            RSSI = RSS(I)
            IF (RSSI.GT.TSS) THEN
               GO TO 60
C
            ELSE
               RSQ(I) = 1.0D0 - RSSI/TSS
               IP = NTERMS(I) + INTP
               IC = N - 2*IP
               IF (IC.LE.0) THEN
                  GO TO 40
C
               ELSE
                  CPI = RSSI/SIGSQ - DBLE(IC)
                  IF (CPI.LT.0.0D0) THEN
                     GO TO 80
C
                  ELSE
                     CP(I) = CPI
                  END IF
               END IF
            END IF
   20    CONTINUE
         GO TO 100
C
   40    NREC = 2
         IERROR = 2
         WRITE (P01REC,FMT=99995) I, N, IP
         GO TO 100
C
   60    IERROR = 3
         NREC = 2
         WRITE (P01REC,FMT=99994) I, RSSI, TSS
         GO TO 100
C
   80    IERROR = 4
         WRITE (P01REC,FMT=99993) I, CPI
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, NMOD.le.0 : NMOD = ',I16)
99998 FORMAT (' ** On entry, SIGSQ.le.0.0 : SIGSQ = ',D13.5)
99997 FORMAT (' ** On entry, TSS.le.0.0 : TSS = ',D13.5)
99996 FORMAT (' ** On entry, MEAN is not valid : MEAN = ',A1)
99995 FORMAT (' ** On entry, number of parameters for model ',I16,' is',
     *       ' too large for N',/'             N = ',I16,' number of p',
     *       'arameters = ',I16)
99994 FORMAT (' ** On entry, value of RSS.gt.TSS :',/'              RS',
     *       'S(',I16,') = ',D13.5,' TSS = ',D13.5)
99993 FORMAT (' ** Value of CP.lt.0.0 : CP(',I16,') = ',D13.5)
      END
