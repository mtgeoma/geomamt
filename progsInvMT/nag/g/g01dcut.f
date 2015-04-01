      INTEGER FUNCTION G01DCU(II,JJ)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATES THE POSITION IN AN ARRAY OF THE (II,JJ)'TH ELEMENT OF
C     A SQUARE MATRIX STORED BY COLUMNS
C
C     ARGUMENTS :
C                 II - ROW POSITION IN ORIGINAL MATRIX
C                 JJ - COLUMN POSITION IN ORIGINAL MATRIX
C
C     .. Scalar Arguments ..
      INTEGER                 II, JJ
C     .. Local Scalars ..
      INTEGER                 I, J
C     .. Executable Statements ..
      I = II
      J = JJ
      IF (I.LE.J) GO TO 20
      I = JJ
      J = II
   20 CONTINUE
      G01DCU = J*(J-1)/2 + I
      RETURN
      END
