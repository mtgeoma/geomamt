      SUBROUTINE G01ASF(PRT,M,N,X,LDX,NSTEPX,NSTEPY,PLOT,LDP,WORK,IWORK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01ASF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDP, LDX, M, NSTEPX, NSTEPY
      CHARACTER         PRT
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(5*M), X(LDX,M)
      INTEGER           IWORK(LDX), N(M)
      CHARACTER         PLOT(LDP,NSTEPX)
C     .. Local Scalars ..
      DOUBLE PRECISION  EXPOS, QPOS, RANGE, WMAX, WMIN, YUNIT
      INTEGER           I, IERR, IEXPOS, IFAIL2, IXUNIT, J, LBODI,
     *                  LENEXT, LLQY, LMEDY, LQLEN, LQPOS, LUQY, MAXY,
     *                  MID, MINY, MTEST, NENDX, NMAX, NOUT
      LOGICAL           IPRINT
      CHARACTER*133     RECORD
C     .. Local Arrays ..
      DOUBLE PRECISION  SCALE(6)
      CHARACTER*9       LABS(6)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01ALF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, MOD, NINT, DBLE
C     .. Executable Statements ..
C
C     Check the input values.
C
      IERR = 0
      IF (NSTEPX.LT.19) THEN
         IERR = 2
         WRITE (REC,FMT=99998) NSTEPX
      ELSE IF (NSTEPY.LT.9) THEN
         IERR = 3
         WRITE (REC,FMT=99997) NSTEPY
      ELSE IF (LDP.LT.NSTEPY) THEN
         IERR = 4
         WRITE (REC,FMT=99996) LDP, NSTEPY
      ELSE
         MTEST = NINT(TWO*((DBLE(NSTEPX)-9.0D0)/(3.0D0*DBLE(M))))
         IF (MTEST.LT.3) THEN
            IERR = 2
            WRITE (REC,FMT=99995) NSTEPX, M
         END IF
      END IF
C
C     Check for valid entry of PRT (PRINT P or NOPRINT N)
C
      IF (PRT.EQ.'P' .OR. PRT.EQ.'p') THEN
         IPRINT = .TRUE.
      ELSE IF (PRT.EQ.'N' .OR. PRT.EQ.'n') THEN
         IPRINT = .FALSE.
      ELSE
         IERR = 5
         WRITE (REC,FMT=99994) PRT
         GO TO 420
      END IF
C
C     Check for LDX being sufficiently large, and N(I) too small.
C
      NMAX = N(1)
      DO 20 I = 1, M
         NMAX = MAX(NMAX,N(I))
   20 CONTINUE
      IF (LDX.LT.NMAX) THEN
         IERR = 6
         WRITE (REC,FMT=99993) LDX, NMAX
      ELSE IF (NMAX.LT.5) THEN
         IERR = 7
         WRITE (REC,FMT=99992)
      END IF
      IF (IERR.EQ.0) THEN
C
C        Initialize PLOT array.
C
         DO 60 I = 1, NSTEPY
            DO 40 J = 1, NSTEPX
               PLOT(I,J) = ' '
   40       CONTINUE
            PLOT(I,9) = ':'
   60    CONTINUE
C
         DO 80 I = 1, M
C
C           For each batch calculate the number summary.
C           If N.lt.5 then omit the batch.
C
            IF (N(I).GE.5) THEN
               IFAIL2 = 0
               CALL G01ALF(N(I),X(1,I),IWORK,WORK((I-1)*5+1),IFAIL2)
C
C              G01ALF can't fail.
C
            ELSE
               IERR = 1
               WRITE (REC,FMT=99999)
            END IF
   80    CONTINUE
C
C        Work out the maximum, and minimum extreme values and hence
C        determine a scale and the appropriate positions.
C
         DO 100 I = 1, M
            IF (N(I).GE.5) THEN
               WMAX = WORK((I-1)*5+5)
               WMIN = WORK((I-1)*5+1)
               GO TO 120
            END IF
  100    CONTINUE
  120    DO 140 I = 1, M
            IF (N(I).GE.5) THEN
               IF (WORK((I-1)*5+1).LT.WMIN) WMIN = WORK((I-1)*5+1)
               IF (WORK((I-1)*5+5).GT.WMAX) WMAX = WORK((I-1)*5+5)
            END IF
  140    CONTINUE
         RANGE = WMAX - WMIN
         IF (RANGE.LE.ZERO) THEN
            IERR = 8
            WRITE (REC,FMT=99991)
            GO TO 420
         END IF
         YUNIT = (NSTEPY-1)/RANGE
         IXUNIT = (NSTEPX-9)/M
C
C        IEXPOS = x-coordinate of left edge of extreme bar.
C        LENEXT = length in array units of the extreme bar.
C        LQPOS = x-coordinate of left edge of quartile bar.
C        LQLEN = length in array units of the quartile bar.
C        MID = x-coordinate of the whiskers.
C
         QPOS = DBLE(IXUNIT)/6.0D0
         EXPOS = TWO*QPOS
         LQPOS = NINT(QPOS)
         LQLEN = NINT(4.0D0*QPOS)
         IF (MOD(LQLEN,2).EQ.0) THEN
            LQLEN = LQLEN - 1
         END IF
         IEXPOS = NINT(EXPOS)
         IF (DBLE(LQPOS).GT.QPOS .AND. DBLE(IEXPOS).LT.EXPOS) THEN
            IEXPOS = IEXPOS + 1
         ELSE IF (DBLE(LQPOS).LT.QPOS .AND. DBLE(IEXPOS).GT.EXPOS) THEN
            IEXPOS = IEXPOS - 1
         END IF
         LENEXT = NINT(HALF*DBLE(LQLEN))
         IF (MOD(LENEXT,2).EQ.0) THEN
            LENEXT = LENEXT - 1
         END IF
         MID = IEXPOS + (LENEXT/2)
C
C        Calculate the appropriate Y-coordinates and plot.
C
         DO 300 I = 1, M
            IF (N(I).GE.5) THEN
               MAXY = NINT(((WMAX-WORK((I-1)*5+5))*YUNIT)+ONE)
               MINY = NINT(((WMAX-WORK((I-1)*5+1))*YUNIT)+ONE)
               LUQY = NINT(((WMAX-WORK((I-1)*5+4))*YUNIT)+ONE)
               LMEDY = NINT(((WMAX-WORK((I-1)*5+3))*YUNIT)+ONE)
               LLQY = NINT(((WMAX-WORK((I-1)*5+2))*YUNIT)+ONE)
C
C              If the height of the plot is less than 5 we can't do it!
C
               IF ((MINY-MAXY).LT.4) THEN
                  DO 160 J = LQPOS + 9, LQPOS + LQLEN + 8
                     PLOT(LMEDY,J) = '*'
  160             CONTINUE
               ELSE
C
C                 Plot extreme bars.
C
                  DO 180 J = IEXPOS + 9, IEXPOS + LENEXT + 8
                     PLOT(MAXY,J) = '-'
                     PLOT(MINY,J) = '-'
  180             CONTINUE
C
C                 Plot the upper whisker.
C
                  DO 200 J = MAXY + 1, LUQY - 1
                     PLOT(J,MID+9) = ':'
  200             CONTINUE
C
C                 Plot the quartile and median bars.
C
                  DO 220 J = LQPOS + 9, LQPOS + LQLEN + 8
                     PLOT(LUQY,J) = '-'
                     PLOT(LMEDY,J) = '-'
                     PLOT(LLQY,J) = '-'
  220             CONTINUE
C
C                 Plot the shell of the box.
C
                  LBODI = LLQY - LUQY
                  DO 240 J = LUQY + 1, LUQY + LBODI - 1
                     PLOT(J,LQPOS+9) = ':'
                     PLOT(J,LQPOS+LQLEN+8) = ':'
  240             CONTINUE
C
C                 Plot the lower whisker.
C
                  DO 260 J = LLQY + 1, MINY - 1
                     PLOT(J,MID+9) = ':'
  260             CONTINUE
               END IF
            ELSE
C
C              Can't plot batch
C
               DO 280 J = -2, 2
                  PLOT(NINT(DBLE(NSTEPY)/TWO)+J,MID+9) = 'X'
  280          CONTINUE
            END IF
C
C           Increment the X-coordinates to move to the next batch
C
            IEXPOS = IEXPOS + IXUNIT
            LQPOS = LQPOS + IXUNIT
            MID = MID + IXUNIT
  300    CONTINUE
C
C        Work out the scale.
C
         SCALE(1) = WMAX
         DO 320 I = 2, 6
            SCALE(I) = SCALE(I-1) - (RANGE/5.0D0)
  320    CONTINUE
         MAXY = WMAX
         WRITE (LABS(1),FMT=99990) SCALE(1)
         WRITE (LABS(6),FMT=99990) SCALE(6)
         DO 340 J = 1, 9
            PLOT(1,J) = LABS(1) (J:J)
            PLOT(NSTEPY,J) = LABS(6) (J:J)
  340    CONTINUE
         DO 380 I = 5, 2, -1
            WRITE (LABS(I),FMT=99990) SCALE(I)
            DO 360 J = 1, 9
               PLOT(NINT((I-1)*(DBLE(NSTEPY)/5.0D0)),J) = LABS(I) (J:J)
  360       CONTINUE
  380    CONTINUE
C
C        Now write out the character arrays if requested.
C
         IF (IPRINT) THEN
            CALL X04ABF(0,NOUT)
            NENDX = MIN(NSTEPX,132)
            DO 400 I = 1, NSTEPY
               WRITE (RECORD,FMT=99989) (PLOT(I,J),J=1,NENDX)
               CALL X04BAF(NOUT,RECORD)
  400       CONTINUE
         END IF
      END IF
  420 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
C
99999 FORMAT (' ** On entry, at least one element of N .lt. 5')
99998 FORMAT (' ** On entry, NSTEPX .lt. 19 :  NSTEPX = ',I13)
99997 FORMAT (' ** On entry, NSTEPY .lt. 9 :  NSTEPY = ',I13)
99996 FORMAT (' ** On entry, LDP.lt.NSTEPY : LDP = ',I16,' NSTEPY = ',
     *       I16)
99995 FORMAT (' ** On entry, NSTEPX = ',I16,' is too small for',I16,
     *       ' plots.')
99994 FORMAT (' ** On entry, PRT is an invalid character :  PRT = ',A1)
99993 FORMAT (' ** On entry, LDX = ',I16,' :  LDX must be .ge. ',I16)
99992 FORMAT (' ** On entry, all batches have N .lt. 5.')
99991 FORMAT (' ** On entry, the data values are all identical.')
99990 FORMAT (D8.1,'+')
99989 FORMAT (1X,132A1)
      END
