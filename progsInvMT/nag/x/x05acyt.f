      LOGICAL FUNCTION X05ACY(STR1,STR2)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Compares two strings for equality, ignoring case. If the strings
C     are of unequal length, the excess characters in the longer string
C     are ignored. Thus, 'April' compares equal to 'APRIL FOOL'.
C
C     X05ACY assumes that the value ICHAR (X) - ICHAR(x) is equal
C     for all upper case alphabetic characters X and their lower
C     case equivalents x. This is true for both ASCII and EBCDIC.
C     The alphabetic characters are not assumed to be stored
C     contiguously.
C
C     .. Scalar Arguments ..
      CHARACTER*(*)           STR1, STR2
C     .. Local Scalars ..
      INTEGER                 DIFF, I, IC, LACODE, LL, UACODE, UZCODE
      CHARACTER*80            A, B
C     .. Intrinsic Functions ..
      INTRINSIC               CHAR, ICHAR, LEN, MIN
C     .. Executable Statements ..
      LL = MIN(LEN(STR1),LEN(STR2))
      A = STR1(1:LL)
      B = STR2(1:LL)
      UACODE = ICHAR('A')
      UZCODE = ICHAR('Z')
      LACODE = ICHAR('a')
      DIFF = UACODE - LACODE
C
      DO 20 I = 1, LL
         IC = ICHAR(A(I:I))
         IF (IC.GE.UACODE .AND. IC.LE.UZCODE) THEN
            A(I:I) = CHAR(IC-DIFF)
         END IF
   20 CONTINUE
C
      DO 40 I = 1, LL
         IC = ICHAR(B(I:I))
         IF (IC.GE.UACODE .AND. IC.LE.UZCODE) THEN
            B(I:I) = CHAR(IC-DIFF)
         END IF
   40 CONTINUE
C
      X05ACY = A .EQ. B
      RETURN
      END
