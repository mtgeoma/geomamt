      INTEGER FUNCTION F07ZAY(NAME)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1993.
C
C     F07ZAY returns a unique positive integer code
C     corresponding to a six-letter NAG routine name
C     given in NAME. If NAME is not recognised, 0 is
C     returned.
C
C     Modified at Mark 16 to allow calls from F08 routines.
C
C     .. Scalar Arguments ..
      CHARACTER*6             NAME
C     .. Local Scalars ..
      INTEGER                 J, K
      CHARACTER               NAME4, NAME5
      CHARACTER*3             NAME13
C     .. Executable Statements ..
C
      IF (NAME(3:3).EQ.'7' .OR. NAME(3:3).EQ.'8') THEN
         NAME13 = NAME(1:3)
         NAME4 = NAME(4:4)
         NAME5 = NAME(5:5)
      ELSE
         NAME13 = NAME(4:6)
         NAME4 = NAME(1:1)
         NAME5 = NAME(2:2)
      END IF
C
      IF (NAME13.EQ.'F07') THEN
C
         IF (NAME4.EQ.'A') THEN
            J = 0
         ELSE IF (NAME4.EQ.'B') THEN
            J = 1
         ELSE IF (NAME4.EQ.'F') THEN
            J = 2
         ELSE IF (NAME4.EQ.'H') THEN
            J = 3
         ELSE IF (NAME4.EQ.'M') THEN
            J = 4
         ELSE IF (NAME4.EQ.'N') THEN
            J = 5
         ELSE IF (NAME4.EQ.'T') THEN
            J = 6
         ELSE
            J = -1
         END IF
C
         IF (NAME5.EQ.'D') THEN
            K = 0
         ELSE IF (NAME5.EQ.'J') THEN
            K = 1
         ELSE IF (NAME5.EQ.'R') THEN
            K = 2
         ELSE IF (NAME5.EQ.'W') THEN
            K = 3
         ELSE
            K = -1
         END IF
C
         IF (J.LT.0 .OR. K.LT.0) THEN
            F07ZAY = 0
         ELSE
C           F07ZAY is in the range 1-28 for F07 routines.
            F07ZAY = 1 + 4*J + K
         END IF
C
      ELSE IF (NAME13.EQ.'F08') THEN
C
         IF (NAME4.EQ.'A') THEN
            J = 0
         ELSE IF (NAME4.EQ.'C') THEN
            J = 1
         ELSE IF (NAME4.EQ.'F') THEN
            J = 2
         ELSE IF (NAME4.EQ.'J') THEN
            J = 3
         ELSE IF (NAME4.EQ.'K') THEN
            J = 4
         ELSE IF (NAME4.EQ.'N') THEN
            J = 5
         ELSE IF (NAME4.EQ.'P') THEN
            J = 6
         ELSE IF (NAME4.EQ.'S') THEN
            J = 7
         ELSE
            J = -1
         END IF
C
         IF (NAME5.EQ.'E') THEN
            K = 0
         ELSE IF (NAME5.EQ.'F') THEN
            K = 1
         ELSE IF (NAME5.EQ.'G') THEN
            K = 2
         ELSE IF (NAME5.EQ.'H') THEN
            K = 3
         ELSE IF (NAME5.EQ.'J') THEN
            K = 4
         ELSE IF (NAME5.EQ.'K') THEN
            K = 5
         ELSE IF (NAME5.EQ.'S') THEN
            K = 6
         ELSE IF (NAME5.EQ.'T') THEN
            K = 7
         ELSE IF (NAME5.EQ.'U') THEN
            K = 8
         ELSE IF (NAME5.EQ.'V') THEN
            K = 9
         ELSE IF (NAME5.EQ.'W') THEN
            K = 10
         ELSE IF (NAME5.EQ.'X') THEN
            K = 11
         ELSE
            K = -1
         END IF
C
         IF (J.LT.0 .OR. K.LT.0) THEN
            F07ZAY = 0
         ELSE
C           F07ZAY is in the range 29-124 for F08 routines.
            F07ZAY = 29 + 12*J + K
         END IF
C
      ELSE
         F07ZAY = 0
      END IF
C
      RETURN
C
      END
