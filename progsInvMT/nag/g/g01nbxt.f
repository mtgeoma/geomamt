      DOUBLE PRECISION FUNCTION G01NBX(LCODE,IS,NK,TR1,TR2,PP,TRACE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C      Based on routine CALCRA by Magnus and Pesaran
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, TWO, FOUR, ZERO
      PARAMETER                        (ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0,
     *                                 ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 TRACE
      INTEGER                          IS, LCODE
C     .. Array Arguments ..
      DOUBLE PRECISION                 PP(*), TR1(*), TR2(*)
      INTEGER                          NK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 PR, SUM, SUM1
      INTEGER                          I, IJ, J, JJ, L, NI, NJ, NL
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE
C     .. Executable Statements ..
      PR = ONE
      DO 20 J = 1, IS
         NJ = NK(J)
         IF (NJ.NE.0) PR = PR*(TR1(J)**NJ)
   20 CONTINUE
      IF (LCODE.EQ.0) THEN
         G01NBX = ZERO
      ELSE IF (LCODE.EQ.1) THEN
         G01NBX = PR
      ELSE
         G01NBX = PR
         SUM = ZERO
         SUM1 = ZERO
         JJ = 0
         IJ = 0
         DO 120 J = 1, IS
            JJ = JJ + J
            NJ = NK(J)
            PR = ONE
            DO 40 L = 1, IS
               NL = NK(L)
               IF (L.EQ.J) NL = NL - 1
               IF (NL.GT.0) PR = PR*(TR1(L)**NL)
   40       CONTINUE
            IF (NJ.NE.0) SUM = SUM + DBLE(J*NJ)*PR*TR2(J)
            IF (LCODE.NE.2) THEN
               PR = ONE
               DO 60 L = 1, IS
                  NL = NK(L)
                  IF (L.EQ.J) NL = NL - 2
                  IF (NL.GT.0) PR = PR*(TR1(L)**NL)
   60          CONTINUE
               IF (NJ.GT.1) SUM1 = SUM1 + DBLE(J*J*NJ*(NJ-1))*PR*PP(JJ)
               DO 100 I = 1, J - 1
                  IJ = IJ + 1
                  NI = NK(I)
                  PR = ONE
                  DO 80 L = 1, IS
                     NL = NK(L)
                     IF (L.EQ.I .OR. L.EQ.J) NL = NL - 1
                     IF (NL.GT.0) PR = PR*(TR1(L)**NL)
   80             CONTINUE
                  IF (NI.GT.0 .AND. NJ.GT.0) SUM1 = SUM1 +
     *                TWO*DBLE(I*J*NI*NJ)*PR*PP(IJ)
  100          CONTINUE
               IJ = IJ + 1
            END IF
  120    CONTINUE
         IF (LCODE.EQ.2) THEN
            G01NBX = G01NBX*TRACE + TWO*SUM
         ELSE
            G01NBX = G01NBX*TRACE + TWO*SUM + FOUR*SUM1
         END IF
      END IF
      RETURN
      END
