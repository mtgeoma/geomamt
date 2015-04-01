      SUBROUTINE D02NUF(NEQ,NEQMAX,JCEVAL,NWKJAC,IA,NIA,JA,NJA,JACPVT,
     *                  NJCPVT,SENS,U,ETA,LBLOCK,ISPLIT,RWORK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-528 (FEB 1987).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NUF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ETA, SENS, U
      INTEGER           IFAIL, ISPLIT, NEQ, NEQMAX, NIA, NJA, NJCPVT,
     *                  NWKJAC
      LOGICAL           LBLOCK
      CHARACTER*1       JCEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(50+4*NEQMAX)
      INTEGER           IA(NIA), JA(NJA), JACPVT(NJCPVT)
C     .. Scalars in Common ..
      CHARACTER*6       LAOPTN
C     .. Local Scalars ..
      INTEGER           I, IAPI, IAPIM1, IDEV, IERR, IPIAN, IPJAN, J, K,
     *                  KCOL, KMAX, KMIN, L, LESTJC, LESTWK, LNJCRQ,
     *                  LNJREQ, LNWKJC, NMOSS
      LOGICAL           LFULL, REPORT
      CHARACTER*1       JCEVL1
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04UDU, X04AAF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Common blocks ..
      COMMON            /AD02NC/LAOPTN
C     .. Save statement ..
      SAVE              /AD02NC/
C     .. Executable Statements ..
      LAOPTN = 'SPARSE'
      IERR = 0
      REPORT = IFAIL .LE. 0
      CALL X04AAF(0,IDEV)
      JCEVL1 = JCEVAL
C     ENSURE UPPER CASE
      CALL E04UDU(JCEVL1)
      IF (JCEVL1.EQ.'S' .OR. JCEVL1.EQ.'D') THEN
         NMOSS = 0
      ELSE IF (JCEVL1.EQ.'A') THEN
         NMOSS = 1
      ELSE IF (JCEVL1.EQ.'N') THEN
         NMOSS = 2
      ELSE IF (JCEVL1.EQ.'F') THEN
         NMOSS = 3
      ELSE
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99999) JCEVL1
            CALL X04BAF(IDEV,REC(1))
            WRITE (REC(2),FMT=99998)
            CALL X04BAF(IDEV,REC(2))
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         RWORK(48) = DBLE(NMOSS)
         RWORK(34) = 3.0D0
         RWORK(35) = 1.0D0
      END IF
      IF (NEQ.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99997) NEQ
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (NEQMAX.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99996) NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (NEQ.GT.NEQMAX) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99995) NEQ, NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (SENS.NE.0.0D0) THEN
         RWORK(43) = ABS(SENS)
      ELSE
         RWORK(43) = 100.0D0*X02AJF()
      END IF
      IF (ISPLIT.LT.1 .OR. ISPLIT.GT.99) THEN
         RWORK(36) = 8.0D0/11.0D0
         RWORK(37) = 0.0D0
      ELSE
         RWORK(36) = DBLE(ISPLIT)/100.0D0
         RWORK(37) = 1.0D0
      END IF
      IF (U.LT.0.0D0 .OR. U.GT.0.9999D0) THEN
         RWORK(38) = 0.1D0
      ELSE
         RWORK(38) = U
      END IF
      IF (ETA.GT.1.0D0) THEN
         RWORK(39) = 2.0D0
      ELSE IF (ETA.LE.0.0D0) THEN
         RWORK(39) = 1.0D-4
      ELSE
         RWORK(39) = ETA
      END IF
      RWORK(41) = 1.0D0
      IF (LBLOCK) THEN
         RWORK(42) = 1.0D0
      ELSE
         RWORK(42) = 0.0D0
      END IF
C
      IF (IERR.EQ.0) THEN
         IF (NMOSS.EQ.1 .OR. NMOSS.EQ.2) THEN
            CALL X04ABF(0,IDEV)
            LESTJC = 20*NEQMAX
            IF (LESTJC.GT.NJCPVT) THEN
               WRITE (REC(1),FMT=99981) NJCPVT, LESTJC
               CALL X04BAF(IDEV,REC(1))
            END IF
            LESTWK = 4*NEQMAX
            IF (LESTWK.GT.NWKJAC) THEN
               WRITE (REC(1),FMT=99980) NWKJAC, LESTWK
               CALL X04BAF(IDEV,REC(1))
            END IF
            GO TO 100
         END IF
      END IF
C
      IF (NIA.LT.NEQ+1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99994) NIA, NEQ + 1
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      LNJREQ = IA(NEQ+1) - 1
      IF (LNJREQ.LT.NEQ) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99993) LNJREQ, NEQ
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (NJA.LT.LNJREQ) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99992) NJA, LNJREQ
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      LNJCRQ = 3*LNJREQ + 14*NEQ + 1
      IF (LNJCRQ.GT.NJCPVT) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99991) NJCPVT, LNJCRQ
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      LNWKJC = LNJREQ + 2*NEQ
      IF (LNWKJC.GT.NWKJAC) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99982) NWKJAC, LNWKJC
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (IA(1).NE.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99990) IA(1)
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      JACPVT(1) = 1
      DO 20 I = 2, NEQ + 1
         IAPI = IA(I)
         IAPIM1 = IA(I-1)
         IF (IAPI.LE.IAPIM1) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(1),FMT=99989) I, IAPI
               CALL X04BAF(IDEV,REC(1))
               WRITE (REC(2),FMT=99988) I - 1, IAPIM1
               CALL X04BAF(IDEV,REC(2))
            END IF
         ELSE IF (IAPI.GT.IAPIM1+NEQ) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(1),FMT=99987) I, IAPI, I - 1
               CALL X04BAF(IDEV,REC(1))
               WRITE (REC(2),FMT=99986) IAPIM1, NEQ
               CALL X04BAF(IDEV,REC(2))
            END IF
         ELSE
            JACPVT(I) = IA(I)
         END IF
   20 CONTINUE
      KCOL = NJCPVT - NEQ
      IPIAN = 1
      IPJAN = IPIAN + NEQ + 1
      KMIN = 1
      DO 80 I = 1, NEQ
         KMAX = IA(I+1) - 1
         DO 40 J = 1, NEQ
            JACPVT(KCOL+J) = 0
   40    CONTINUE
         DO 60 K = KMIN, KMAX
            L = JA(K)
            IF (L.LT.1 .OR. L.GT.NEQ) THEN
               IERR = 1
               IF (REPORT) THEN
                  WRITE (REC(1),FMT=99985) K, L
                  CALL X04BAF(IDEV,REC(1))
                  WRITE (REC(2),FMT=99984) NEQ
                  CALL X04BAF(IDEV,REC(2))
               END IF
               GO TO 60
            END IF
            LFULL = JACPVT(KCOL+L) .EQ. 1
            IF (LFULL) THEN
               IERR = 1
               IF (REPORT) THEN
                  WRITE (REC(3),FMT=99983) I, L
                  CALL X04BAF(IDEV,REC(3))
               END IF
            ELSE
               JACPVT(KCOL+L) = 1
               JACPVT(K+IPJAN-1) = L
            END IF
   60    CONTINUE
         KMIN = KMAX + 1
   80 CONTINUE
  100 CONTINUE
      IF (IERR.EQ.0) THEN
         RWORK(46) = DBLE(NWKJAC)
         RWORK(47) = DBLE(NJCPVT)
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
99999 FORMAT (' ** D02NUF - JCEVAL(1:1)(=',A1,') IS NOT')
99998 FORMAT ('           ''N'' OR ''S'' OR ''A'' OR ''F'' OR ''D'' **')
99997 FORMAT (' ** D02NUF - NEQ(=',I16,') .LT. 1 ** ')
99996 FORMAT (' ** D02NUF - NEQMAX(=',I16,') .LT. 1 ** ')
99995 FORMAT (' ** D02NUF - NEQ(=',I16,') .GT. NEQMAX(=',I16,') ** ')
99994 FORMAT (' ** D02NUF - NIA(=',I16,') .LT. NEQ+1(=',I16,') ** ')
99993 FORMAT (' ** D02NUF - IA(NEQ+1)-1(=',I16,') .LT. NEQ(=',I16,
     *  ' ** ')
99992 FORMAT (' ** D02NUF - NJA(=',I16,').LT.IA(NEQ+1)-1(=',I16,') ** ')
99991 FORMAT (' ** D02NUF - NJCPVT(=',I16,').LT.LENGTH REQUIRED(=',I16,
     *  ') ** ')
99990 FORMAT (' ** D02NUF - IA(1)(=',I16,') IS NOT EQUAL TO 1 ** ')
99989 FORMAT (' ** D02NUF - IA(',I16,')(=',I16,') .LT.')
99988 FORMAT ('             IA(',I16,')(=',I16,') ** ')
99987 FORMAT (' ** D02NUF - IA(',I16,')(=',I16,') - IA(',I16,')')
99986 FORMAT ('             (=',I16,') .GT. NEQ(=',I16,') ** ')
99985 FORMAT (' ** D02NUF - JA(',I16,')(=',I16,') .LT. 1  OR ')
99984 FORMAT ('             .GT. NEQ(=',I16,') ** ')
99983 FORMAT (' ** D02NUF - DUPLICATE ENTRY IN ROW ',I16,',COL',I16,
     *  ' ** ')
99982 FORMAT (' ** D02NUF - NWKJAC(=',I16,').LT.LENGTH REQUIRED(=',I16,
     *  ') ** ')
99981 FORMAT (' ** D02NUF - NJCPVT(=',I16,') .LT. GUESS(=',I16,') ** ')
99980 FORMAT (' ** D02NUF - NWKJAC(=',I16,') .LT. GUESS(=',I16,') ** ')
      END
