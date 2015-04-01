      SUBROUTINE G01ARV(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,IWORK,
     *                  FROM,TO,HIEND,SLBRK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Print the LO or HI stem for a Stem-and-Leaf display.
C     The logical variable  HIEND  is  .TRUE.  if we are to print
C     the  HI  stem, and  .FALSE.  if the  LO  stem is to be printed.
C     IWORK()  contains  N  sorted and scaled data values.  Each has the
C     form SS...SL, where the one's digit is the leaf.
C     FROM, TO  are pointer into  IWORK() delimiting the values to be
C     placed on the  HI  or  LO  stem.
C     SLBRK  is the character position on the page of the blank column
C     between stems and leaves.
C
C     .. Scalar Arguments ..
      INTEGER           FROM, LDP, MAXPTR, OUTPTR, PMAX, PSUB, SLBRK, TO
      LOGICAL           HIEND
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      CHARACTER         CHARS(*), PLOT(LDP,*)
C     .. Local Scalars ..
      INTEGER           CHBL, CHCOMA, CHH, CHI, CHL, CHO, I, LHMAX,
     *                  NWID, OPOS
C     .. External Functions ..
      INTEGER           G01ARR
      EXTERNAL          G01ARR
C     .. External Subroutines ..
      EXTERNAL          G01ARY, G01ARZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Data statements ..
      DATA              CHH, CHI, CHL, CHO, CHCOMA, CHBL/13, 14, 15, 17,
     *                  28, 22/
C     .. Executable Statements ..
C
      OPOS = SLBRK - 3
      IF (HIEND) THEN
         OUTPTR = 1
         MAXPTR = 1
         PSUB = PSUB + 2
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,CHH)
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CHI)
      ELSE
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,CHL)
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CHO)
      END IF
      CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,SLBRK,CHBL)
      NWID = MAX(G01ARR(IWORK(FROM)),G01ARR(IWORK(TO)))
      LHMAX = PMAX - NWID - 2
      DO 20 I = FROM, TO
         CALL G01ARY(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,IWORK(I),
     *               NWID)
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CHCOMA)
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,0,CHBL)
         IF (OUTPTR.GT.LHMAX) THEN
            OUTPTR = 1
            MAXPTR = 1
            PSUB = PSUB + 1
            CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,SLBRK,
     *                  CHBL)
         END IF
   20 CONTINUE
C
C        But dont print the final comma
C
      OPOS = MAXPTR - 1
      CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,OPOS,CHBL)
      OUTPTR = 1
      MAXPTR = 1
      IF ( .NOT. HIEND) PSUB = PSUB + 1
   40 RETURN
      END
