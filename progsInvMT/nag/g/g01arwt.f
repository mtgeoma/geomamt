      SUBROUTINE G01ARW(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,EPSI,
     *                  MAXINT,W,IWORK,N,PT1,PT2,CUT,IADJH,HI,RANK,
     *                  MEDYET,DEPWID,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Compute and print the depth for the current line.
C
C     W()       := the  N  sorted data values.
C     IWORK()   := the scaled version of W().
C     PT1, PT2  := pointers into IWORK() and W().
C     ON ENTRY, PT1 = PT2 point to the first datum not yet printed.
C     ON EXIT,  PT1 points to 1st datum on next line, PT2 is unchanged.
C     CUT       := the largest value on the current (positive) line,
C               or the smallest value above the current (negative) line.
C     IADJH     points to the high adjacent value in W() and IWORK()
C     HI        := the greatest value being displayed
C     RANK      := a running total of the rank from the low end.
C     ON EXIT, RANK  is updated to include count for the current line.
C     MEDYET    := a logical flag, set  .TRUE.  when the median
C                  value has been processed.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, TWO
      PARAMETER         (ZERO=0.0D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSI, MAXINT
      INTEGER           CUT, DEPWID, HI, IADJH, IFAIL, LDP, MAXPTR, N,
     *                  OUTPTR, PMAX, PSUB, PT1, PT2, RANK
      LOGICAL           MEDYET
C     .. Array Arguments ..
      DOUBLE PRECISION  W(*)
      INTEGER           IWORK(*)
      CHARACTER         CHARS(*), PLOT(LDP,*)
C     .. Local Scalars ..
      INTEGER           CHLPAR, CHRPAR, DEPTH, LEFCNT, NWID, OPOS, PTX,
     *                  PTZ
C     .. External Functions ..
      INTEGER           G01ARR, G01ARS
      EXTERNAL          G01ARR, G01ARS
C     .. External Subroutines ..
      EXTERNAL          G01ARY, G01ARZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              CHLPAR, CHRPAR/26, 27/
C     .. Executable Statements ..
C
      PTX = PT1
      DO 20 PT1 = PTX, IADJH
         IF (IWORK(PT1).GT.CUT) GO TO 40
         IF ((CUT.GE.0) .AND. (IWORK(PT1).EQ.CUT)) GO TO 40
   20 CONTINUE
      GO TO 140
   40 IF (CUT.NE.0) THEN
         GO TO 160
C
C        Zero Cut:  if data all .le. 0, all zeroes go on "-0" stem
C
      ELSE IF (HI.GT.0) THEN
C
C        Both +0 and -0 stems -- share the zeroes between them
C
C        First check for numbers rounded to zero--true -0s
         DO 60 PTZ = PT1, N
            IF (W(PTZ).GE.ZERO) GO TO 80
   60    CONTINUE
   80    PT1 = PTZ
         DO 100 PTZ = PT1, N
            IF (W(PTZ).GT.ZERO) GO TO 120
  100    CONTINUE
  120    PT1 = PT1 + G01ARS(EPSI,MAXINT,DBLE(PTZ-PT1)/TWO,IFAIL)
         GO TO 160
      END IF
C
C     Last data value: if we fall past--point past it for consistency.
C
  140 PT1 = IADJH + 1
C
C     Compute and print depth
C
  160 LEFCNT = PT1 - PT2
      RANK = RANK + LEFCNT
C
C     Case: where is the median?
C
      IF (MEDYET) THEN
C
C        Case 1: past the median
C
         DEPTH = N - (RANK-LEFCNT)
      ELSE
         IF (DBLE(RANK).EQ.DBLE(N)/TWO) THEN
C
C           Case 2: median falls between stems at this point
C
            MEDYET = .TRUE.
         ELSE IF (DBLE(RANK).GE.DBLE(N+1)/TWO) THEN
C
C           Case 3: median is on the current line
C
            NWID = G01ARR(LEFCNT)
            OPOS = DEPWID - NWID
            CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,
     *                  CHLPAR)
            CALL G01ARY(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,LEFCNT,
     *                  NWID)
            CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CHRPAR)
            MEDYET = .TRUE.
            GO TO 180
         END IF
C
C        Case 4: not up to median yet
C
         DEPTH = RANK
      END IF
C
C     Print the depth, if it hasn't been done yet
C
      NWID = G01ARR(DEPTH)
      OPOS = DEPWID - NWID + 2
      CALL G01ARY(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,DEPTH,
     *            NWID)
  180 RETURN
      END
