      SUBROUTINE C06FAZ(PTS,PMAX,TWOGRP,TFACT,RFACT,IERROR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COPY REORDERING FACTORING PROGRAMME
C     .. Scalar Arguments ..
      INTEGER           IERROR, PMAX, PTS, TWOGRP
C     .. Array Arguments ..
      INTEGER           RFACT(21), TFACT(21)
C     .. Local Scalars ..
      INTEGER           F, J, JJ, N, NEST, P, PTWO, Q
C     .. Local Arrays ..
      INTEGER           PP(10), QQ(20)
C     .. Data statements ..
      DATA              NEST/20/
C     .. Executable Statements ..
      N = PTS
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
      GO TO 20
C
  100 CONTINUE
      IF (P.LT.1) GO TO 140
      DO 120 J = 1, P
         JJ = P + 1 - J
         TFACT(J) = PP(JJ)
         JJ = P + Q + J
         TFACT(JJ) = PP(J)
  120 CONTINUE
  140 CONTINUE
      IF (Q.LT.1) GO TO 180
      DO 160 J = 1, Q
         JJ = P + J
         TFACT(JJ) = QQ(J)
  160 CONTINUE
  180 CONTINUE
      JJ = 2*P + Q
      TFACT(JJ+1) = 0
      RFACT(JJ+1) = 0
      DO 200 J = 1, JJ
         RFACT(J) = TFACT(J)
  200 CONTINUE
      IF (JJ.EQ.1) RFACT(1) = 0
      PTWO = 1
      J = 0
  220 CONTINUE
      J = J + 1
      IF (TFACT(J).EQ.0) GO TO 260
      IF (TFACT(J).NE.2) GO TO 220
      PTWO = PTWO*2
      TFACT(J) = 1
      IF (PTWO.GE.TWOGRP) GO TO 240
      IF (TFACT(J+1).EQ.2) GO TO 220
  240 CONTINUE
      TFACT(J) = PTWO
      PTWO = 1
      GO TO 220
  260 CONTINUE
      IERROR = 0
      GO TO 320
  280 IERROR = 1
      GO TO 320
C
  300 IERROR = 2
  320 RETURN
      END
