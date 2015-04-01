      SUBROUTINE G01ART(EPSI,MAXINT,HI,LO,MAXP,PTOTL,FRACT,UNIT,INUNIT,
     *                  NPW,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Find a simple data-units value to assign to one plot position
C     in one dimension of a plot.  A plot position is typically
C     one character position horizontally, or one line vertically.
C
C     ON ENTRY:
C     HI, LO  := high and low edges of the data range to be plotted.
C     MAXP    := the maximum number of plot positions allowed in this
C                dimension of the plot.
C
C     ON EXIT:
C     PTOTL      holds the total number of plot positions to be used in
C                this dimension.  (must be  .le.  MAXP.)
C     FRACT   := the mantissa of the nice position width.
C                it is selected from the numbers in NICNOS.
C     UNIT    := an integer power of 10 such that NPW = FRACT * UNIT.
C     NPW     := the nice position width.  One plot position width
C                will represent  a data-space distance of  NPW.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, TEN
      INTEGER           NN
      PARAMETER         (ZERO=0.0D0,TEN=10.0D0,NN=4)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSI, FRACT, HI, LO, MAXINT, NPW, UNIT
      INTEGER           IFAIL, MAXP, PTOTL
      LOGICAL           INUNIT
C     .. Local Scalars ..
      DOUBLE PRECISION  APRXW, Z
      INTEGER           FLOOR, I
C     .. Local Arrays ..
      DOUBLE PRECISION  NICNOS(NN)
C     .. External Functions ..
      INTEGER           G01ARS
      EXTERNAL          G01ARS
C     .. Intrinsic Functions ..
      INTRINSIC         INT, LOG10, DBLE
C     .. Data statements ..
      DATA              NICNOS/1.0D0, 2.0D0, 5.0D0, 10.0D0/
C     .. Executable Statements ..
C
C     FLOOR is the largest integer not exceeding Z.
      IF (MAXP.GT.0) THEN
         APRXW = (HI-LO)/DBLE(MAXP)
         IF ( .NOT. INUNIT) THEN
            Z = LOG10(APRXW)
            FLOOR = INT(Z)
            IF (Z.LT.ZERO .AND. Z.NE.DBLE(FLOOR)) FLOOR = FLOOR - 1
            UNIT = TEN**FLOOR
         END IF
         FRACT = APRXW/UNIT
C
         DO 20 I = 1, NN
            IF (FRACT.LE.NICNOS(I)) GO TO 40
   20    CONTINUE
C
   40    CONTINUE
         IF (I.GT.NN) I = NN
         FRACT = NICNOS(I)
         NPW = FRACT*UNIT
         PTOTL = G01ARS(EPSI,MAXINT,HI/NPW,IFAIL) - G01ARS(EPSI,MAXINT,
     *           LO/NPW,IFAIL) + 1
C
         IF (IFAIL.EQ.0) THEN
C
C           If minus zero position possible and SGN(HI) .ne. SGN(LO),
C           allow it.
C
            IF (HI*LO.LT.ZERO .OR. HI.EQ.ZERO) PTOTL = PTOTL + 1
C
C           PTOTL positions required with this width -- few enough?
C
            IF (PTOTL.GT.MAXP) THEN
C
C              Too many positions needed, so increase NPW.
C
               IF (I.EQ.NN .AND. INUNIT) THEN
C                 Number of lines required is greater than MAXP
                  IFAIL = 3
                  GO TO 60
               ELSE IF (I.EQ.NN) THEN
                  I = 1
                  UNIT = UNIT*TEN
               ELSE
                  I = I + 1
               END IF
               GO TO 40
            END IF
         END IF
      ELSE
         IFAIL = 3
      END IF
   60 RETURN
      END
