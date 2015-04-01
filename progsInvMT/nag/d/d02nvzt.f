      SUBROUTINE D02NVZ(NEQMAX,NY2DIM,MAXORD,METH,PETZLD,CONST,TCRIT,
     *                  HMIN,HMAX,H0,MAXSTP,MXHNIL,NORM,RWORK,IERR,
     *                  REPORT)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-529 (FEB 1987).
C     MARK 13B REVISED. IER-655 (AUG 1988).
C
C     SPECIAL SETUP ROUTINE FOR SHELL
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H0, HMAX, HMIN, TCRIT
      INTEGER           IERR, MAXORD, MAXSTP, METH, MXHNIL, NEQMAX,
     *                  NY2DIM
      LOGICAL           PETZLD, REPORT
      CHARACTER*1       NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  CONST(6), RWORK(50+4*NEQMAX)
C     .. Local Scalars ..
      INTEGER           I, IDEV, NMETH
      CHARACTER*1       NORM1
C     .. Local Arrays ..
      DOUBLE PRECISION  CONSTX(6)
      INTEGER           MORD(2)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          E04UDU, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              CONSTX/2.0D0, 1.0D1, 1.0D3, 1.2D0, 1.3D0, 1.4D0/
      DATA              MORD/12, 5/
C     .. Executable Statements ..
      CALL X04AAF(0,IDEV)
      IF (NEQMAX.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99992) NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         RWORK(1) = TCRIT
         RWORK(2) = DBLE(MAXSTP)
         IF (MXHNIL.GT.0) THEN
            RWORK(3) = DBLE(MXHNIL)
         ELSE
            RWORK(3) = 10.0D0
         END IF
         RWORK(5) = H0
         RWORK(6) = HMAX
         RWORK(7) = HMIN
         RWORK(8) = 1.0D0
         RWORK(40) = 0.0D0
      END IF
      NORM1 = NORM
C     ENSURE UPPER CASE
      CALL E04UDU(NORM1)
      IF (NORM1.EQ.'M') THEN
         RWORK(49) = 2.0D0
      ELSE IF (NORM1.EQ.'A' .OR. NORM1.EQ.'D') THEN
         RWORK(49) = 1.0D0
      ELSE
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99999) NORM1
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      RWORK(21) = 1.0D0
      RWORK(22) = 1.0D0
      DO 20 I = 1, 6
         IF (CONST(I).LT.0.0D0) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(2),FMT=99998) I, CONST(I)
               CALL X04BAF(IDEV,REC(2))
            END IF
         ELSE IF (CONST(I).EQ.0.0D0) THEN
            CONST(I) = CONSTX(I)
         END IF
         IF (CONST(I).LT.1.0D0 .AND. I.GT.1) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(2),FMT=99997) I, CONST(I)
               CALL X04BAF(IDEV,REC(2))
            END IF
         ELSE
            RWORK(22+I) = CONST(I)
         END IF
   20 CONTINUE
      IF (CONST(1).GT.CONST(2) .OR. CONST(2).GT.CONST(3)) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(2),FMT=99996)
            CALL X04BAF(IDEV,REC(2))
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         IF (METH.EQ.10 .OR. METH.EQ.20) THEN
            NMETH = METH/10
            RWORK(29) = 1.0D0
            RWORK(30) = 0.0D0
         ELSE
            IF (METH.NE.2) METH = 1
            NMETH = METH
            IF (PETZLD) THEN
               RWORK(30) = 1.0D0
            ELSE
               RWORK(30) = 0.0D0
            END IF
            RWORK(29) = 0.0D0
         END IF
         RWORK(31) = DBLE(NMETH)
         IF (MAXORD.LE.0 .OR. MAXORD.GT.MORD(NMETH)) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(3),FMT=99995) MAXORD, MORD(NMETH)
               CALL X04BAF(IDEV,REC(3))
            END IF
         END IF
      END IF
      RWORK(32) = DBLE(MAXORD)
      IF (NY2DIM.LT.MAXORD+1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(3),FMT=99994) MAXORD
            CALL X04BAF(IDEV,REC(3))
            WRITE (REC(4),FMT=99993) MAXORD + 1, NY2DIM
            CALL X04BAF(IDEV,REC(4))
         END IF
      END IF
      IF (IERR.EQ.0) RWORK(33) = DBLE(NY2DIM)
      RETURN
C
99999 FORMAT (' ** D02NVF - NORM(1:1)(=',A1,') IS NOT ''A'' OR ''M'' O',
     *  'R ''D'' **')
99998 FORMAT (' ** D02NVF - CONST(',I1,') = ',1P,D12.5,' .LT. 0.0 **')
99997 FORMAT (' ** D02NVF - CONST(',I1,') = ',1P,D12.5,' .LT. 1.0 **')
99996 FORMAT (' ** D02NVF - CONST(1).GT.CONST(2) OR CONST(2).GT.CONST(',
     *  '3) **')
99995 FORMAT (' ** D02NVF - MAXORD = ',I16,' DOES NOT LIE IN (1,',I2,
     *  ') **')
99994 FORMAT (' ** D02NVF - MAXORD = ',I2,
     *  ', HENCE NY2DIM SHOULD EQUAL,')
99993 FORMAT ('    AT LEAST ',I2,' BUT USER SUPPLIED NY2DIM =',I16,
     *  ' **')
99992 FORMAT (' ** D02NVF - NEQMAX(=',I16,') .LT. 1 **')
      END
