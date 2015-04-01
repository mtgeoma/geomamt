      INTEGER FUNCTION G01ARR(I)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Find the number of characters needed to print I
C
C     .. Scalar Arguments ..
      INTEGER                 I
C     .. Local Scalars ..
      INTEGER                 IA, IQ, ND
C     .. Intrinsic Functions ..
      INTRINSIC               ABS
C     .. Executable Statements ..
C
      IA = ABS(I)
C
      ND = 1
      IF (I.LT.0) ND = 2
   20 CONTINUE
      IQ = IA/10
      IF (IQ.NE.0) THEN
         IA = IQ
         ND = ND + 1
         GO TO 20
      END IF
C
      G01ARR = ND
C
      END
