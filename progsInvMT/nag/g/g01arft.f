      SUBROUTINE G01ARF(RANGE,PRT,N,Y,NSTEPX,NSTEPY,UNIT,PLOT,LDP,LINES,
     *                  SORTY,IWORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Produce a Stem-and-Leaf display of the data in Y()
C
C     IWORK() is an integer work array. SORTY() is a real work array.
C
C     UNIT       is the leaf unit used.  Is converted to a power of ten.
C     NSTEPX     is the maximum width of the plot horizontally
C     NSTEPY     is the maximum number of lines to feature the display.
C     LDP        is the leading dimension of the plotting array.
C     EPSI       is the machine-related epsilon.
C     IMAXIN     is the maximum permitted integer value
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01ARF')
      DOUBLE PRECISION  ZERO, ONE, TEN
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TEN=10.0D0)
      INTEGER           CMAX
      PARAMETER         (CMAX=29)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  UNIT
      INTEGER           IFAIL, LDP, LINES, N, NSTEPX, NSTEPY
      CHARACTER         PRT, RANGE
C     .. Array Arguments ..
      DOUBLE PRECISION  SORTY(N), Y(N)
      INTEGER           IWORK(N)
      CHARACTER         PLOT(LDP,NSTEPX)
C     .. Local Scalars ..
      DOUBLE PRECISION  ADJH, ADJL, ALOW, EPSI, FRACT, HH, HL, MAXINT,
     *                  MED, MNUNIT, NPW, STEP, X, Z
      INTEGER           ADDLIN, CHSTAR, CUT, DEPWID, FLOOR, HI, I,
     *                  IADJH, IADJL, IERR, IMAXIN, J, K, L, LEAF,
     *                  LINWID, LOW, MAXPTR, NDIGIT, NENDX, NHI, NLINHI,
     *                  NLINLO, NLINS, NLMAX, NLO, NOUT, NVLINH, NVLINL,
     *                  NWIDH, NWIDL, OUTPTR, PLTWID, PT1, PT2, PVMAX,
     *                  RANGI, RANK, SLBRK, SPACNT, STEM, STMWID, XH, XL
      LOGICAL           INUNIT, IPRINT, MEDYET, NEGNOW, XTREMS
      CHARACTER         BLANK
      CHARACTER*133     RECORD
C     .. Local Arrays ..
      CHARACTER         CHARS(CMAX)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           G01ARR, G01ARS, P01ABF, X02BBF
      EXTERNAL          X02AJF, G01ARR, G01ARS, P01ABF, X02BBF
C     .. External Subroutines ..
      EXTERNAL          G01ARQ, G01ART, G01ARU, G01ARV, G01ARW, G01ARX,
     *                  G01ARY, G01ARZ, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, LOG10, MAX, MIN, MOD, DBLE
C     .. Data statements ..
C     DATA definitions: useful characters and the scaling options
      DATA              CHSTAR, BLANK/25, ' '/
      DATA              CHARS/'0', '1', '2', '3', '4', '5', '6', '7',
     *                  '8', '9', 'E', 'F', 'H', 'I', 'L', 'N', 'O',
     *                  'P', 'R', 'S', 'T', ' ', '+', '-', '*', '(',
     *                  ')', ',', '.'/
C     .. Executable Statements ..
C
      IERR = 1
      IF (N.LT.2) THEN
         WRITE (REC,FMT=99999) N
      ELSE IF (NSTEPX.LT.35) THEN
         WRITE (REC,FMT=99998) NSTEPX
      ELSE IF (NSTEPY.GT.0 .AND. NSTEPY.LT.5) THEN
         WRITE (REC,FMT=99997) NSTEPY
      ELSE IF (LDP.LT.5) THEN
         WRITE (REC,FMT=99996) LDP
      ELSE IF (LDP.LT.NSTEPY) THEN
         WRITE (REC,FMT=99995) LDP, NSTEPY
      ELSE
         IERR = 0
      END IF
      IF (PRT.EQ.'P' .OR. PRT.EQ.'p') THEN
         IPRINT = .TRUE.
      ELSE IF (PRT.EQ.'N' .OR. PRT.EQ.'n') THEN
         IPRINT = .FALSE.
      ELSE
         IERR = 2
         WRITE (REC,FMT=99994) PRT
      END IF
      IF (RANGE.EQ.'E' .OR. RANGE.EQ.'e') THEN
         XTREMS = .TRUE.
      ELSE IF (RANGE.EQ.'F' .OR. RANGE.EQ.'f') THEN
         XTREMS = .FALSE.
      ELSE
         IERR = 2
         WRITE (REC,FMT=99993) RANGE
      END IF
      IF (IERR.EQ.0) THEN
         OUTPTR = 1
         MAXPTR = 1
         INUNIT = .FALSE.
         EPSI = X02AJF()
         IMAXIN = X02BBF(X)
         MAXINT = DBLE(IMAXIN)
         NDIGIT = ABS(LOG10(MAXINT)) + 1
C
         DO 40 J = 1, NSTEPX
            DO 20 I = 1, LDP
               PLOT(I,J) = BLANK
   20       CONTINUE
   40    CONTINUE
C
C        Sort Y in SORTY and get summary information
C
         DO 60 I = 1, N
            SORTY(I) = Y(I)
   60    CONTINUE
         CALL G01ARQ(SORTY,N,MED,HL,HH,ADJL,ADJH,IADJL,IADJH,STEP)
C
C        Find nice line width for plot
C        If adjacent values equal or user demands it, fake the adjacent
C        values to be the extremes
C
         IF (ADJH.LE.ADJL .OR. XTREMS) THEN
            IADJL = 1
            IADJH = N
            ADJL = SORTY(IADJL)
            ADJH = SORTY(IADJH)
         END IF
C
         DEPWID = G01ARR(N)
         IF (UNIT.GT.ZERO) THEN
            INUNIT = .TRUE.
            Z = LOG10(UNIT)
            FLOOR = INT(Z)
            IF (Z.LT.ZERO .AND. Z.NE.DBLE(FLOOR)) FLOOR = FLOOR - 1
            UNIT = TEN**FLOOR
         END IF
C
         IF (ADJH.LE.ADJL) THEN
C           Even if all values are equal we can produce a display
            NLINS = 1
            FRACT = TEN
            IF ( .NOT. INUNIT) THEN
               ALOW = ABS(ADJL)
               IF (ALOW.EQ.ZERO) THEN
                  UNIT = ONE
               ELSE
                  K = LOG10(ALOW)
                  L = LOG10(EPSI)
                  IF (K.GE.0) K = K + 1
                  UNIT = TEN**(K+L+1)
                  MNUNIT = ALOW/MAXINT
                  IF (UNIT.LT.MNUNIT) THEN
                     Z = LOG10(MNUNIT) + 1
                     FLOOR = INT(Z)
                     IF (Z.LT.ZERO .AND. Z.NE.DBLE(FLOOR))
     *                   FLOOR = FLOOR - 1
                     UNIT = TEN**FLOOR
                  END IF
               END IF
            END IF
            IWORK(1) = G01ARS(EPSI,MAXINT,ADJL/UNIT,IERR)
C
            DO 80 I = 2, N
               IWORK(I) = IWORK(1)
   80       CONTINUE
            GO TO 160
         END IF
C
C        Work out number of lines for HI - LO
C
         NLINLO = 0
         NLINHI = 0
         NLO = IADJL - 1
         NHI = N - IADJH
         IF (NLO.GT.0 .OR. NHI.GT.0) THEN
            IF (INUNIT) THEN
               XL = G01ARS(EPSI,MAXINT,SORTY(1)/UNIT,IERR)
               NWIDL = G01ARR(XL)
               XH = G01ARS(EPSI,MAXINT,SORTY(N)/UNIT,IERR)
               NWIDH = G01ARR(XH)
               XL = G01ARS(EPSI,MAXINT,ADJL/UNIT,IERR)
               XH = G01ARS(EPSI,MAXINT,ADJH/UNIT,IERR)
               STMWID = MAX(G01ARR(XL),G01ARR(XH)) - 1
            ELSE
               NWIDL = NDIGIT
               NWIDH = NDIGIT
               XL = G01ARS(EPSI,MAXINT,ADJL,IERR)
               XH = G01ARS(EPSI,MAXINT,ADJH,IERR)
               STMWID = MAX(G01ARR(XL),G01ARR(XH))
            END IF
            RANGI = NSTEPX - DEPWID - STMWID - 3
            IF (NLO.GT.0) THEN
               NVLINL = RANGI/(NWIDL+2)
               IF (NVLINL.EQ.0) THEN
                  NLINLO = NLO + 1
               ELSE
                  NLINLO = NLO/NVLINL + 2
               END IF
            END IF
            IF (NHI.GT.0) THEN
               NVLINH = RANGI/(NWIDH+2)
               IF (NVLINH.EQ.0) THEN
                  NLINHI = NHI + 2
               ELSE
                  NLINHI = NHI/NVLINH + 3
               END IF
            END IF
         END IF
         ADDLIN = NLINLO + NLINHI
C        Calculate NLMAX
         IF (NSTEPY.GT.0) THEN
            NLMAX = NSTEPY - ADDLIN - 4
            PVMAX = NSTEPY
         ELSE
            NLMAX = G01ARS(EPSI,MAXINT,TEN*LOG10(DBLE(IADJH-IADJL+1)),
     *              IERR)
            IF (LDP.LT.(NLMAX+ADDLIN+4)) NLMAX = LDP - ADDLIN - 4
            PVMAX = LDP
         END IF
         IF (IERR.EQ.4) THEN
            WRITE (REC,FMT=99992)
            GO TO 240
         END IF
C        Calculate UNIT and FRACT
         CALL G01ART(EPSI,MAXINT,ADJH,ADJL,NLMAX,NLINS,FRACT,UNIT,
     *               INUNIT,NPW,IERR)
         IF (IERR.EQ.3) THEN
            WRITE (REC,FMT=99991) (NLINS+ADDLIN+4), PVMAX
            GO TO 240
         ELSE IF (IERR.EQ.4) THEN
            WRITE (REC,FMT=99992)
            GO TO 240
         END IF
C
C        Rescale everything according to UNIT.  Hereafter everything is
C        integer, and data are of the form   SS...SL(.)
C        NOTE: G01ARS performs epsilon adjustments for correct rounding,
C        and checks the real number is not too large for an integer
C        variable.
C
         DO 100 I = 1, N
            IWORK(I) = G01ARS(EPSI,MAXINT,SORTY(I)/UNIT,IERR)
  100    CONTINUE
         IF (IERR.EQ.4) THEN
            WRITE (REC,FMT=99992)
            GO TO 240
         END IF
C
         IF ( .NOT. INUNIT) THEN
C           If all leaves are 0, should be in one-line-per-stem format
            DO 120 I = IADJL, IADJH
               IF (MOD(IWORK(I),10).NE.0) GO TO 160
  120       CONTINUE
            FRACT = ONE
            UNIT = UNIT*TEN
            NPW = FRACT*UNIT
            NLINS = G01ARS(EPSI,MAXINT,ADJH/NPW,IERR) - G01ARS(EPSI,
     *              MAXINT,ADJL/NPW,IERR) + 1
            IF (ADJH*ADJL.LT.ZERO .OR. ADJH.EQ.ZERO) NLINS = NLINS + 1
            DO 140 I = 1, N
               IWORK(I) = IWORK(I)/10
  140       CONTINUE
         END IF
  160    LOW = IWORK(IADJL)
         HI = IWORK(IADJH)
C
C        Linewidth now is NICEWIDTH/UNIT = FRACT
C
         LINWID = G01ARS(EPSI,MAXINT,FRACT,IERR)
C
C        Setup :- find width of plotting region,
C                 Stem-Leaf break position, etc.
C
         STMWID = MAX(G01ARR(HI),G01ARR(LOW))
         SLBRK = STMWID + DEPWID + 5
         IF (LINWID.NE.10) SLBRK = SLBRK + 2
         PLTWID = NSTEPX - SLBRK - 2
C
         CALL G01ARU(PLOT,LDP,NSTEPX,OUTPTR,MAXPTR,CHARS,EPSI,MAXINT,
     *               UNIT,IERR)
         IF (IERR.EQ.4) THEN
            WRITE (REC,FMT=99992)
            GO TO 240
         END IF
C
C        Print values below low adjacent value on "LO" stem.
C
         RANK = IADJL - 1
         IF (RANK.EQ.0) THEN
            LINES = 4
         ELSE
            LINES = 5
            OUTPTR = 1
            MAXPTR = 1
            CALL G01ARV(PLOT,LDP,LINES,NSTEPX,OUTPTR,MAXPTR,CHARS,IWORK,
     *                  1,RANK,.FALSE.,SLBRK)
         END IF
C
C        Initialize for main part of display.
C        Initial settings are to line before first one printed
C        FLOOR is the largest integer not exceeding Z.
C
         Z = (ONE+EPSI)*DBLE(LOW)/DBLE(LINWID)
         FLOOR = INT(Z)
         IF (Z.LT.ZERO .AND. Z.NE.DBLE(FLOOR)) FLOOR = FLOOR - 1
         CUT = FLOOR*LINWID
         NEGNOW = .TRUE.
         STEM = CUT
         IF (LOW.GE.0) THEN
C           First stem positive
            NEGNOW = .FALSE.
            STEM = CUT - LINWID
         END IF
         MEDYET = .FALSE.
C
C        2 pointers are used.  PT1 counts first for depths, PT2 follows
C        for leaf printing.  Both are initialized one point early.
C
         PT1 = IADJL
         PT2 = PT1
C
C        Main loop.  For each line
C
         DO 200 J = 1, NLINS
C
C           Variable used:
C                CUT = first number on next line of positive stems, but
C                    = last number on current line of negative stems
C               STEM = inner (near zero) edge of current line
C             SPACNT   counts spaces used on this line
C
C           Step to next line
            OUTPTR = 1
            MAXPTR = 1
            LINES = LINES + 1
C
            CUT = CUT + LINWID
            IF (STEM.EQ.0 .AND. NEGNOW) THEN
               NEGNOW = .FALSE.
            ELSE
               STEM = STEM + LINWID
            END IF
C
C           Newline -- initialize count of spaces used
C
            SPACNT = 0
C
C           Find and print depth
C
            CALL G01ARW(PLOT,LDP,LINES,NSTEPX,OUTPTR,MAXPTR,CHARS,EPSI,
     *                  MAXINT,SORTY,IWORK,N,PT1,PT2,CUT,IADJH,HI,RANK,
     *                  MEDYET,DEPWID,IERR)
            IF (IERR.EQ.4) THEN
               WRITE (REC,FMT=99992)
               GO TO 240
            END IF
C
C           Print stem label
C
            CALL G01ARX(PLOT,LDP,LINES,NSTEPX,OUTPTR,MAXPTR,CHARS,STEM,
     *                  LINWID,NEGNOW,SLBRK)
C
C           Find and print leaves
C
            IF (PT1.NE.PT2) THEN
  180          LEAF = MOD(ABS(IWORK(PT2)-(STEM/10)*10),10)
               CALL G01ARY(PLOT,LDP,LINES,NSTEPX,OUTPTR,MAXPTR,CHARS,0,
     *                     LEAF,1)
               SPACNT = SPACNT + 1
               IF (SPACNT.GE.PLTWID) THEN
C
C                 Line overflows past right edge.  Mark with *
C
                  CALL G01ARZ(PLOT,LDP,LINES,NSTEPX,OUTPTR,MAXPTR,CHARS,
     *                        0,CHSTAR)
                  PT2 = PT1
               ELSE
                  PT2 = PT2 + 1
                  IF (PT2.LT.PT1) GO TO 180
               END IF
C
C              End line
C
            END IF
  200    CONTINUE
C
C        Print values above HI adjacent value on "HI" stem.
C
         IF (PT1.LE.N) THEN
            CALL G01ARV(PLOT,LDP,LINES,NSTEPX,OUTPTR,MAXPTR,CHARS,IWORK,
     *                  PT1,N,.TRUE.,SLBRK)
         END IF
C
         IF (IPRINT) THEN
            CALL X04ABF(0,NOUT)
            NENDX = MIN(132,NSTEPX)
            DO 220 I = 1, LINES
               WRITE (RECORD,FMT=99990) (PLOT(I,J),J=1,NENDX)
               CALL X04BAF(NOUT,RECORD)
  220       CONTINUE
         END IF
      END IF
  240 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (' ** On entry,  N .lt. 2 : N = ',I16)
99998 FORMAT (' ** On entry,  NSTEPX .lt. 35 : NSTEPX = ',I16)
99997 FORMAT (' ** On entry,  NSTEPY .gt. 0  and  .lt. 5 : NSTEPY =',
     *       I16)
99996 FORMAT (' ** On entry,  LDP .lt. 5 : LDP = ',I16)
99995 FORMAT (' ** On entry, LDP .lt. NSTEPY : LDP =',I16,'  NSTEPY =',
     *       I16)
99994 FORMAT (' ** On entry,  PRT is an invalid character: PRT = ',A1)
99993 FORMAT (' ** On entry,  RANGE is an invalid character: RANGE = ',
     *       A1)
99992 FORMAT (' ** A value exceeds maximum allowed for an integer')
99991 FORMAT (' ** Lines needed for display (',I16,') exceed NSTEPY (',
     *       I16,')')
99990 FORMAT (1X,132A1)
      END
