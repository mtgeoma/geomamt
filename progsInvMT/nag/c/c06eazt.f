      SUBROUTINE C06EAZ(PTS,PMAX,TWOGRP,FACTOR,SYM,PSYM,UNSYM,IERROR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     SYMMETRIZED REORDERING FACTORING PROGRAMME
C     .. Scalar Arguments ..
      INTEGER           IERROR, PMAX, PSYM, PTS, TWOGRP
C     .. Array Arguments ..
      INTEGER           FACTOR(21), SYM(21), UNSYM(21)
C     .. Local Scalars ..
      INTEGER           F, J, JJ, N, NEST, P, PTWO, Q, R
C     .. Local Arrays ..
      INTEGER           PP(10), QQ(20)
C     .. Data statements ..
      DATA              NEST/20/
C     .. Executable Statements ..
      N = PTS
      PSYM = 1
      F = 2
      P = 0
      Q = 0
   20 CONTINUE
      IF (N.LE.1) GO TO 100
      DO 40 J = F, PMAX
         IF (N.EQ.(N/J)*J) GO TO 60
   40 CONTINUE
      GO TO 280
   60 CONTINUE
      IF (2*P+Q.GE.NEST) GO TO 300
      F = J
      N = N/F
      IF (N.EQ.(N/F)*F) GO TO 80
      Q = Q + 1
      QQ(Q) = F
      GO TO 20
   80 CONTINUE
      N = N/F
      P = P + 1
      PP(P) = F
      PSYM = PSYM*F
      GO TO 20
C
  100 CONTINUE
      R = 1
      IF (Q.EQ.0) R = 0
      IF (P.LT.1) GO TO 140
      DO 120 J = 1, P
         JJ = P + 1 - J
         SYM(J) = PP(JJ)
         FACTOR(J) = PP(JJ)
         JJ = P + Q + J
         FACTOR(JJ) = PP(J)
         JJ = P + R + J
         SYM(JJ) = PP(J)
  120 CONTINUE
  140 CONTINUE
      IF (Q.LT.1) GO TO 180
      DO 160 J = 1, Q
         JJ = P + J
         UNSYM(J) = QQ(J)
         FACTOR(JJ) = QQ(J)
  160 CONTINUE
      SYM(P+1) = PTS/PSYM**2
  180 CONTINUE
      JJ = 2*P + Q
      FACTOR(JJ+1) = 0
      PTWO = 1
      J = 0
  200 CONTINUE
      J = J + 1
      IF (FACTOR(J).EQ.0) GO TO 240
      IF (FACTOR(J).NE.2) GO TO 200
      PTWO = PTWO*2
      FACTOR(J) = 1
      IF (PTWO.GE.TWOGRP) GO TO 220
      IF (FACTOR(J+1).EQ.2) GO TO 200
  220 CONTINUE
      FACTOR(J) = PTWO
      PTWO = 1
      GO TO 200
  240 CONTINUE
      IF (P.EQ.0) R = 0
      JJ = 2*P + R
      SYM(JJ+1) = 0
      IF (Q.LE.1) Q = 0
      UNSYM(Q+1) = 0
      IERROR = 0
C
  260 CONTINUE
      RETURN
C
  280 CONTINUE
      IERROR = 1
      GO TO 260
C
  300 CONTINUE
      IERROR = 2
      GO TO 260
C
      END
