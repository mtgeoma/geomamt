      SUBROUTINE G01ARU(PLOT,LDP,NSTEPX,OUTPTR,MAXPTR,CHARS,EPSI,MAXINT,
     *                  UNIT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Print the title for a Stem-and-Leaf display.
C
C     PLOT    " STEM-AND-LEAF DISPLAY "
C             " LEAF DIGIT UNIT = " UNIT
C             "  1  2  REPRESENTS " 12*UNIT
C
C     ON ENTRY:
C     UNIT    := the leaf digit unit
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSI, MAXINT, UNIT
      INTEGER           IFAIL, LDP, MAXPTR, NSTEPX, OUTPTR
C     .. Array Arguments ..
      CHARACTER         CHARS(*), PLOT(LDP,*)
C     .. Local Scalars ..
      INTEGER           CH0, CH1, CH2, CHBL, CHE, CHF, CHL, CHMIN, CHN,
     *                  CHP, CHPT, CHR, CHS, CHT, I, IEXPT, J, NUM,
     *                  OWID, PSUB, ZEROS
C     .. Local Arrays ..
      INTEGER           K(3)
      CHARACTER*22      NAME(3)
C     .. External Functions ..
      INTEGER           G01ARR, G01ARS
      EXTERNAL          G01ARR, G01ARS
C     .. External Subroutines ..
      EXTERNAL          G01ARY, G01ARZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG10
C     .. Data statements ..
C
C     Assign the headings to the PLOT array
C
      DATA              CH0, CH1, CH2, CHE, CHF, CHL, CHN, CHP, CHR,
     *                  CHS, CHT, CHBL, CHMIN, CHPT/1, 2, 3, 11, 12, 15,
     *                  16, 18, 19, 20, 21, 22, 24, 29/
      DATA              NAME/'Stem-and-leaf display',
     *                  'Leaf digit unit = ', '1  2  represents  '/
      DATA              K/21, 18, 18/
C     .. Executable Statements ..
      DO 40 I = 1, 3
         DO 20 J = 1, K(I)
            PLOT(I,J) = NAME(I) (J:J)
   20    CONTINUE
   40 CONTINUE
C
      IEXPT = G01ARS(EPSI,MAXINT,LOG10(UNIT),IFAIL)
C
C     If      IEXPT .ge.  0          UNIT .ge. 1.0
C     If      IEXPT .eq. -1          UNIT .eq. 0.1
C     Else    IEXPT .lt.  0          UNIT .le. 0.01
C
C     Plot value of UNIT
C
      PSUB = 2
      OUTPTR = 19
      MAXPTR = 19
      IF (IEXPT.GE.0) THEN
         NUM = G01ARS(EPSI,MAXINT,UNIT,IFAIL)
         IF (IFAIL.NE.0) GO TO 100
         OWID = G01ARR(NUM)
         CALL G01ARY(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,NUM,
     *               OWID)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CHPT)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH0)
      ELSE IF (IEXPT.EQ.(-1)) THEN
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH0)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CHPT)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH1)
      ELSE
         ZEROS = ABS(IEXPT) - 1
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH0)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CHPT)
         IF (ZEROS.NE.0) THEN
            DO 60 I = 1, ZEROS
               CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,
     *                     CH0)
   60       CONTINUE
         END IF
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH1)
      END IF
C
C     Plot what ' 1  2 ' represents
C
      PSUB = 3
      OUTPTR = 19
      MAXPTR = 19
      IF (IEXPT.GE.0) THEN
C        NUM = 12 * UNIT
         NUM = G01ARS(EPSI,MAXINT,12.0D0*UNIT,IFAIL)
         IF (IFAIL.NE.0) GO TO 100
         OWID = G01ARR(NUM)
         CALL G01ARY(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,NUM,
     *               OWID)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CHPT)
      ELSE IF (IEXPT.EQ.(-1)) THEN
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH1)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CHPT)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH2)
      ELSE
         ZEROS = ABS(IEXPT) - 2
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH0)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CHPT)
         IF (ZEROS.NE.0) THEN
            DO 80 I = 1, ZEROS
               CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,
     *                     CH0)
   80       CONTINUE
         END IF
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH1)
         CALL G01ARZ(PLOT,LDP,PSUB,NSTEPX,OUTPTR,MAXPTR,CHARS,0,CH2)
      END IF
  100 RETURN
      END
