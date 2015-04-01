      SUBROUTINE D02MWZ(NEQMAX,NY2DIM,THETA,CONST,TCRIT,HMIN,HMAX,H0,
     *                  MAXSTP,MXHNIL,NORM,RWORK,IERR,REPORT)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H0, HMAX, HMIN, TCRIT, THETA
      INTEGER           IERR, MAXSTP, MXHNIL, NEQMAX, NY2DIM
      LOGICAL           REPORT
      CHARACTER*(*)     NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  CONST(3), RWORK(50+4*NEQMAX)
C     .. Local Scalars ..
      INTEGER           I, IDEV, ISAVE
      CHARACTER*1       NORM1
C     .. Local Arrays ..
      DOUBLE PRECISION  CONSTX(3)
      CHARACTER*80      REC(1)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
C
      DATA              CONSTX/2.0D0, 1.0D1, 1.0D3/
C     .. Executable Statements ..
C
      CALL X04AAF(0,IDEV)
      ISAVE = 0
      RWORK(1) = TCRIT
      RWORK(2) = DBLE(MAXSTP)
      RWORK(3) = DBLE(MXHNIL)
      RWORK(5) = H0
      RWORK(6) = HMAX
      RWORK(7) = HMIN
      RWORK(8) = 1.0D0
      RWORK(40) = 0.0D0
      NORM1 = NORM(1:1)
      IF (NORM1.EQ.'M' .OR. NORM1.EQ.'m') THEN
         RWORK(49) = 2.0D0
      ELSE IF (NORM1.EQ.'A' .OR. NORM1.EQ.'a' .OR. NORM1.EQ.'D' .OR.
     *         NORM1.EQ.'d') THEN
         RWORK(49) = 1.0D0
      ELSE
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99999) NORM1
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      RWORK(21) = 3.0D0
      RWORK(22) = 1.0D0
      DO 20 I = 1, 3
         IF (CONST(I).LT.0.0D0) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(1),FMT=99998) I, I, CONST(I)
               CALL X04BAF(IDEV,REC(1))
            END IF
         ELSE IF (CONST(I).EQ.0.0D0) THEN
            CONST(I) = CONSTX(I)
         END IF
         IF (CONST(I).LT.1.0D0 .AND. I.GT.1) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(1),FMT=99997) I, I, CONST(I)
               CALL X04BAF(IDEV,REC(1))
            END IF
         ELSE
            RWORK(22+I) = CONST(I)
         END IF
   20 CONTINUE
      IF (CONST(1).GE.CONST(2) .OR. CONST(2).GE.CONST(3)) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99996)
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      RWORK(29) = 0.0D0
      RWORK(30) = 0.0D0
      RWORK(33) = DBLE(NY2DIM)
      RETURN
C
99999 FORMAT (' ** On entry, NORM is not valid : NORM = ',A1)
99998 FORMAT (' ** On entry, CONST(',I1,').lt.0.0 : CONST(',I1,') =',
     *       D13.5)
99997 FORMAT (' ** On entry, CONST(',I1,').lt.1.0 : CONST(',I1,') =',
     *       D13.5)
99996 FORMAT (' ** On entry, CONST(1).ge.CONST(2) or CONST(2).ge.CO',
     *       'NST(3)')
      END
