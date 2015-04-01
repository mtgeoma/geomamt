      SUBROUTINE D02QDY(N,TN,HLAST,HNEXT,Y,JCEVAL,JACSTR,ML,MU,COUT,CIN,
     *                  COMM,IMON,CHK1,CHK2,ISAVE,X)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HLAST, HNEXT, TN, X
      INTEGER           IMON, ISAVE, ML, MU, N
      CHARACTER*1       JACSTR
      CHARACTER*6       JCEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  CHK1(N), CHK2(N), CIN(7), COMM(5), COUT(16),
     *                  Y(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  RDUM, XSTEPS
      INTEGER           IDUM77, N1, N2, NJE, NQ, NQU, NRE, NST
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, YTEMP, YTEMP1, YTEMP2
      INTEGER           I, NTEMP, NTEMP1
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, DBLE
C     .. Common blocks ..
      COMMON            /BD02NM/RDUM, NQ, NQU, NST, NRE, NJE, N1, N2,
     *                  IDUM77
      COMMON            /DD02QD/XSTEPS
C     .. Save statement ..
      SAVE              /BD02NM/, /DD02QD/
C     .. Executable Statements ..
      IF (IMON.EQ.1 .OR. IMON.EQ.-1) THEN
         COUT(9) = XSTEPS
         COUT(15) = DBLE(NQU)
         COUT(16) = DBLE(NJE)
         CIN(7) = DBLE(NQ)
         CIN(6) = HNEXT
      END IF
      IF (IMON.NE.1) GO TO 100
      COUT(8) = COUT(8) + 1.0D0
      IF (COUT(8).EQ.1.0D0) THEN
         COUT(1) = HLAST
         COUT(2) = HLAST
      ELSE
         COUT(1) = MIN(COUT(1),ABS(HLAST))
         COUT(2) = MAX(COUT(2),ABS(HLAST))
      END IF
      IF (ABS(HLAST).EQ.COUT(14)) COUT(3) = COUT(3) + 1.0D0
      IF (COUT(8).GT.1.0D0) COUT(5) = COUT(4)
      COUT(4) = TN - HLAST
      IF (COUT(8).EQ.1.0D0) CIN(5) = HLAST
      YTEMP = 0.0D0
      DO 20 I = 1, N
         YTEMP = MAX(YTEMP,ABS(Y(I)))
   20 CONTINUE
      COUT(6) = MAX(COUT(6),YTEMP)
      COUT(7) = MIN(COUT(7),YTEMP)
      IF (COMM(3).EQ.0.0D0) GO TO 60
      DO 40 I = 1, N
         YTEMP1 = Y(I) - CHK1(I)
         YTEMP2 = CHK2(I) - CHK1(I)
         IF (ABS(YTEMP1).LT.ABS(YTEMP2)) CHK2(I) = Y(I)
         IF (SIGN(1.0D0,YTEMP1).EQ.SIGN(1.0D0,YTEMP2)) GO TO 40
         ISAVE = 0
         IMON = -2
         CIN(1) = 3.0D0
         COMM(3) = 0.0D0
         GO TO 120
   40 CONTINUE
C
   60 CONTINUE
      IF (COMM(2).EQ.0.0D0) GO TO 80
      IF (COMM(2).LE.YTEMP) THEN
         ISAVE = 0
         IMON = -2
         CIN(1) = 4.0D0
         COMM(2) = 0.0D0
         GO TO 120
      END IF
   80 CONTINUE
      IF (COMM(4).EQ.0.0D0) GO TO 100
      IF (COMM(4).GT.0.0D0) THEN
         ISAVE = 0
         CIN(1) = 5.0D0
         GO TO 120
      ELSE
         IF (SIGN(1.D0,TN-COMM(5)).NE.SIGN(1.D0,COUT(4)-COMM(5))) THEN
            ISAVE = 0
            IMON = -2
            CIN(1) = 6.0D0
            COMM(4) = 0.0D0
            GO TO 120
         END IF
      END IF
  100 CONTINUE
      IF (COMM(1).GT.0.0D0) THEN
         NTEMP = NRE
         IF (JCEVAL(1:1).EQ.'A') THEN
            IF (JACSTR.EQ.'B') THEN
               NTEMP1 = ML + MU + 1
            ELSE
               NTEMP1 = N
            END IF
            NTEMP = NTEMP + NTEMP1*NJE
         END IF
         TEMP = COMM(1) - DBLE(NTEMP)
         IF (TEMP.LE.0.0D0) THEN
            IMON = -2
            COMM(1) = MIN(COMM(1),-1.0D0)
            ISAVE = 4
            IF (COUT(8).EQ.0.0D0) ISAVE = 5
            GO TO 120
         END IF
      END IF
      RETURN
  120 CONTINUE
      X = TN
      RETURN
      END
