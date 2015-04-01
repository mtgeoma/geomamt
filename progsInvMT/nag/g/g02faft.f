      SUBROUTINE G02FAF(N,IP,NRES,RES,H,RMS,SRES,LDS,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-887 (NOV 1990).
C
C     CALCULATES STANDARDIZED RESIDUALS FROM RESIDUALS
C     AND DIAGONAL OF HAT MATRIX
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02FAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RMS
      INTEGER           IFAIL, IP, LDS, N, NRES
C     .. Array Arguments ..
      DOUBLE PRECISION  H(NRES), RES(NRES), SRES(LDS,4)
C     .. Local Scalars ..
      DOUBLE PRECISION  P, PN, PN1, R, RD, RS, SPN
      INTEGER           I, IERROR, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (IP.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) IP
      ELSE IF (N-1.LE.IP) THEN
         WRITE (P01REC(1),FMT=99997) IP, N
      ELSE IF (RMS.LE.0.0D0) THEN
         WRITE (P01REC(1),FMT=99996) RMS
      ELSE IF (NRES.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) NRES
      ELSE IF (NRES.GT.N) THEN
         WRITE (P01REC(1),FMT=99992) NRES, N
      ELSE IF (LDS.LT.NRES) THEN
         WRITE (P01REC(1),FMT=99995) LDS, NRES
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
C
C        CALCULATE CONSTANTS
C
         P = DBLE(IP)
         PN = DBLE(N-IP)
         PN1 = PN - 1.0D0
         SPN = SQRT(PN/P)
         DO 20 I = 1, NRES
            IF (RES(I).NE.0.0D0) THEN
               IF (H(I).LE.0.0D0 .OR. H(I).GE.1.0D0) THEN
                  GO TO 40
C
               ELSE
                  R = (1.0D0-H(I))*RMS
                  RD = RES(I)/SQRT(R)
                  R = RES(I)*RES(I)/R
                  RS = (PN-R)
                  IF (RS.LE.0.0D0) THEN
                     GO TO 60
C
                  ELSE
                     RS = RD*SQRT(PN1/RS)
C
C                    SRES(I,1) STANDARDIZED RESIUAL
C
                     SRES(I,1) = RD
C
C                    SRES(I,2) STUDENTIZED RESIDUAL
C
                     SRES(I,2) = RS
C
C                    SRES(I,3) COOKK'S D
C
                     SRES(I,3) = R*H(I)/((1.0D0-H(I))*P)
C
C                    SRES(I,4) ATKINSON'S T
C
                     SRES(I,4) = SQRT(H(I)/(1.0D0-H(I)))*SPN*RS
                  END IF
               END IF
            ELSE
               SRES(I,1) = 0.0D0
               SRES(I,2) = 0.0D0
               SRES(I,3) = 0.0D0
               SRES(I,4) = 0.0D0
            END IF
   20    CONTINUE
         GO TO 80
C
   40    IERROR = 2
         NREC = 2
         WRITE (P01REC,FMT=99994) I, H(I)
         GO TO 80
C
   60    IERROR = 3
         NREC = 2
         WRITE (P01REC,FMT=99993) I, RES(I), RMS
      END IF
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, NRES.lt.1 : NRES = ',I16)
99998 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99997 FORMAT (' ** On entry, N-1.le.IP : IP = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, RMS.le.0 : RMS = ',D13.5)
99995 FORMAT (' ** On entry, LDS.lt.NRES : LDS = ',I16,' NRES = ',I16)
99994 FORMAT (' ** On entry, value in H .le. 0.0 or .ge. 1.0:',/'     ',
     *       '                 H(',I16,') = ',D13.5)
99993 FORMAT (' ** On entry, a value in RES is too large for given RMS',
     *       /'             RES(',I16,') = ',D13.5,' RMS = ',D13.5)
99992 FORMAT (' ** On entry, NRES.gt.N : NRES = ',I16,' N = ',I16)
      END
