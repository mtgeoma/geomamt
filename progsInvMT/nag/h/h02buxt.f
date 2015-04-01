      LOGICAL FUNCTION H02BUX(STRING)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Scalar Arguments ..
      CHARACTER*8             STRING
C     .. Local Scalars ..
      INTEGER                 I
      LOGICAL                 VALID
      CHARACTER               ITEM
C     .. Intrinsic Functions ..
      INTRINSIC               LGE, LLE
C     .. Executable Statements ..
      H02BUX = .TRUE.
C
      IF (STRING(1:1).EQ.' ') THEN
         H02BUX = .FALSE.
         RETURN
      END IF
C
      DO 20 I = 1, 8
         ITEM = STRING(I:I)
         VALID = (LGE(ITEM,'a') .AND. LLE(ITEM,'z')) .OR. (LGE(ITEM,'A')
     *            .AND. LLE(ITEM,'Z')) .OR. (LGE(ITEM,'0')
     *           .AND. LLE(ITEM,'9')) .OR. ITEM .EQ. ' ' .OR. ITEM .EQ.
     *           '$' .OR. ITEM .EQ. '*' .OR. ITEM .EQ. ':' .OR.
     *           ITEM .EQ. '+' .OR. ITEM .EQ. '-' .OR. ITEM .EQ. '.'
         IF ( .NOT. VALID) THEN
            H02BUX = .FALSE.
            RETURN
         END IF
   20 CONTINUE
C
      RETURN
      END
