      CHARACTER*80 FUNCTION X05ACX(STRING,WRDLEN)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     X05ACX returns all the characters of STRING before the first
C     space character, through the function name. The number of
C     characters is returned through WRDLEN. In addition, STRING
C     is stripped of the characters returned in X05ACX.
C     If STRING is entirely blank, X05ACX returns ' ' with WRDLEN = 0.
C
C     .. Scalar Arguments ..
      INTEGER                      WRDLEN
      CHARACTER*(*)                STRING
C     .. Local Scalars ..
      INTEGER                      I, K, LL
C     .. Intrinsic Functions ..
      INTRINSIC                    LEN
C     .. Executable Statements ..
      LL = LEN(STRING)
      IF (STRING.EQ.' ') THEN
         X05ACX = ' '
         WRDLEN = 0
      ELSE
         I = 1
   20    IF (STRING(I:I).EQ.' ') THEN
            I = I + 1
            GO TO 20
         END IF
         K = I + 1
   40    IF (STRING(K:K).NE.' ' .AND. K.LE.LL) THEN
            K = K + 1
            GO TO 40
         END IF
         X05ACX = STRING(I:K-1)
         WRDLEN = K - I
         IF (K.LE.LL) THEN
            STRING = STRING(K:LL)
         ELSE
            STRING = ' '
         END IF
      END IF
      RETURN
      END
