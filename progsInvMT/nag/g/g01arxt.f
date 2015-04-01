      SUBROUTINE G01ARX(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,STEM,
     *                  LINWID,NEGNOW,SLBRK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Compute and print the stem.
C
C     ON ENTRY:
C     STEM    := the inner (near zero) edge of the current line.
C     LINWID  := the number of possible different leaf digits.
C     SLBRK   := the character position on the page of the blank
C                column between stems and leaves.
C     NEGNOW     is  .TRUE. if the current line is negative.
C
C     .. Scalar Arguments ..
      INTEGER           LDP, LINWID, MAXPTR, OUTPTR, PMAX, PSUB, SLBRK,
     *                  STEM
      LOGICAL           NEGNOW
C     .. Array Arguments ..
      CHARACTER         CHARS(*), PLOT(LDP,*)
C     .. Local Scalars ..
      INTEGER           CH0, CHBL, CHMIN, CHPLUS, CHPT, CHSTAR, I,
     *                  LEFDIG, NSTEM, NWID, OCHR, OPOS
C     .. Local Arrays ..
      INTEGER           CH5STM(5)
C     .. External Functions ..
      INTEGER           G01ARR
      EXTERNAL          G01ARR
C     .. External Subroutines ..
      EXTERNAL          G01ARY, G01ARZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              CH0, CHBL, CHPLUS, CHMIN/1, 22, 23, 24/
      DATA              CHSTAR, CHPT/25, 29/
      DATA              CH5STM/25, 21, 12, 20, 29/
C     .. Executable Statements ..
C
      NSTEM = STEM/10
      LEFDIG = ABS(STEM-NSTEM*10)
      NWID = G01ARR(NSTEM)
C
C     Case: how many possible digits/line ( = LINWID)
C
      IF (LINWID.EQ.2) THEN
C
C        Case 1: 2 possible digits/line; 5 lines/stem.
C
         IF (NSTEM.NE.0) THEN
            OPOS = SLBRK - NWID - 2
         ELSE
C           Plus or Minus zero
            OPOS = SLBRK - 4
            IF (NEGNOW) THEN
               CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,
     *                     CHMIN)
            ELSE
               CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,
     *                     CHPLUS)
            END IF
            OPOS = OPOS + 1
         END IF
         CALL G01ARY(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,NSTEM,
     *               NWID)
         I = LEFDIG/2 + 1
         OCHR = CH5STM(I)
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,OCHR)
      ELSE IF (LINWID.EQ.5) THEN
C
C        Case 2: 5 possible digits/line; 2 lines/stem
C
         OPOS = SLBRK - NWID - 1
         IF (NSTEM.EQ.0) THEN
C
C           -0*   print the sign (it appears automatically otherwise)
C
            OPOS = SLBRK - 3
            IF (NEGNOW) THEN
               CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,
     *                     CHMIN)
            ELSE
               CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,
     *                     CHPLUS)
            END IF
         END IF
         OPOS = SLBRK - NWID - 1
         CALL G01ARY(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,NSTEM,
     *               NWID)
         IF (LEFDIG.LT.5) THEN
            CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CHSTAR)
         ELSE
            CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CHPT)
         END IF
C
C        Case 3: 10 possible digits/leaf; 1 line/stem.
C
      ELSE IF ((NSTEM.NE.0) .OR. .NOT. NEGNOW) THEN
         OPOS = SLBRK - NWID - 1
         CALL G01ARY(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,NSTEM,
     *               NWID)
      ELSE
         OPOS = SLBRK - 3
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,CHMIN)
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CH0)
      END IF
      CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,SLBRK,CHBL)
   20 RETURN
      END
