      INTEGER FUNCTION X05AAZ(STR)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Converts the character string STR into an integer.
C     If STR contains any character not in the range '0' .. '9', or
C     if the result would overflow, then X05AAZ returns 0. (For this
C     reason, we do not simply use an internal READ statement to do
C     the conversion).
C
C     X05AAZ assumes that the characters '0', '1', .. '9' are stored
C     contiguously in the machine's character set. This is true for
C     ASCII and EBCDIC.
C
C     .. Scalar Arguments ..
      CHARACTER*(*)           STR
C     .. Local Scalars ..
      INTEGER                 I, K, MAXINT, ZERO
      LOGICAL                 OUTRNG
      CHARACTER               CH
C     .. External Functions ..
      INTEGER                 X02BBF
      EXTERNAL                X02BBF
C     .. Intrinsic Functions ..
      INTRINSIC               ICHAR, LEN
C     .. Executable Statements ..
      ZERO = ICHAR('0')
      MAXINT = X02BBF(1.0D0)
      CH = STR(1:1)
C     Count blanks as zeros.
      IF (CH.EQ.' ') CH = '0'
      K = ICHAR(CH) - ZERO
      IF (K.LT.0 .OR. K.GT.9) THEN
         OUTRNG = .TRUE.
      ELSE
         OUTRNG = .FALSE.
         X05AAZ = K
         DO 20 I = 2, LEN(STR)
            CH = STR(I:I)
            IF (CH.EQ.' ') CH = '0'
            K = ICHAR(CH) - ZERO
            IF (K.LT.0 .OR. K.GT.9) THEN
               OUTRNG = .TRUE.
            ELSE IF (X05AAZ.GT.(MAXINT-K)/10) THEN
               OUTRNG = .TRUE.
            ELSE
               X05AAZ = 10*X05AAZ + K
            END IF
   20    CONTINUE
      END IF
      IF (OUTRNG) X05AAZ = 0
      RETURN
      END
