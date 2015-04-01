      SUBROUTINE C06ECY(X,Y,PTS,SYM,PSYM,UNSYM)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DOUBLE IN PLACE REORDERING PROGRAMME
C     .. Scalar Arguments ..
      INTEGER           PSYM, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
      INTEGER           SYM(21), UNSYM(21)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           DK, I, II, IL, J, JJ, JL, K, KK, KS, LK, MODS,
     *                  MULT, NEST, PUNSYM, TEST
      LOGICAL           ONEMOD
C     .. Local Arrays ..
      INTEGER           MODULO(20)
C     .. External Subroutines ..
      EXTERNAL          C06ECX
C     .. Data statements ..
      DATA              NEST/20/
C     .. Executable Statements ..
      CALL C06ECX(X,Y,PTS,SYM)
C
      IF (UNSYM(1).EQ.0) GO TO 280
      PUNSYM = PTS/PSYM**2
      MULT = PUNSYM/UNSYM(1)
      TEST = (UNSYM(1)*UNSYM(2)-1)*MULT*PSYM
      LK = MULT
      DK = MULT
      DO 20 K = 2, NEST
         IF (UNSYM(K).EQ.0) GO TO 40
         LK = LK*UNSYM(K-1)
         DK = DK/UNSYM(K)
         MODULO(K) = (LK-DK)*PSYM
         MODS = K
   20 CONTINUE
   40 CONTINUE
      ONEMOD = MODS .LT. 3
      IF (ONEMOD) GO TO 80
      K = (MODS+3)/2
      DO 60 J = 3, K
         JJ = MODS + 3 - J
         KK = MODULO(J)
         MODULO(J) = MODULO(JJ)
         MODULO(JJ) = KK
   60 CONTINUE
   80 CONTINUE
      JL = (PUNSYM-3)*PSYM
      KS = PUNSYM*PSYM
C
      DO 260 J = PSYM, JL, PSYM
         JJ = J
C
  100    CONTINUE
         JJ = JJ*MULT
         IF (ONEMOD) GO TO 140
         DO 120 I = 3, MODS
            JJ = JJ - (JJ/MODULO(I))*MODULO(I)
  120    CONTINUE
  140    CONTINUE
         IF (JJ.GE.TEST) GO TO 160
         JJ = JJ - (JJ/MODULO(2))*MODULO(2)
         GO TO 180
  160    CONTINUE
         JJ = JJ - (JJ/MODULO(2))*MODULO(2) + MODULO(2)
  180    CONTINUE
         IF (JJ.LT.J) GO TO 100
C
         IF (JJ.EQ.J) GO TO 240
         LK = JJ - J
         II = J + 1
         IL = J + PSYM
         DO 220 I = II, IL
            DO 200 K = I, PTS, KS
               KK = K + LK
               T = X(K)
               X(K) = X(KK)
               X(KK) = T
               T = Y(K)
               Y(K) = Y(KK)
               Y(KK) = T
  200       CONTINUE
  220    CONTINUE
  240    CONTINUE
  260 CONTINUE
C
  280 CONTINUE
      RETURN
      END
