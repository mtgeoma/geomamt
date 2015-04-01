      SUBROUTINE C06ECX(X,Y,PTS,SYM)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11 REVISED. IER-444 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DOUBLE SYMMETRIC REORDERING PROGRAMME
C     EQUIVALENCE (I1,I(1)), (K1,K(1)), (L1,L(1))
C     EQUIVALENCE (I2,I(2)), (K2,K(2)), (L2,L(2))
C     EQUIVALENCE (I3,I(3)), (K3,K(3)), (L3,L(3))
C     EQUIVALENCE (I4,I(4)), (K4,K(4)), (L4,L(4))
C     EQUIVALENCE (I5,I(5)), (K5,K(5)), (L5,L(5))
C     EQUIVALENCE (I6,I(6)), (K6,K(6)), (L6,L(6))
C     EQUIVALENCE (I7,I(7)), (K7,K(7)), (L7,L(7))
C     EQUIVALENCE (I8,I(8)), (K8,K(8)), (L8,L(8))
C     EQUIVALENCE (I9,I(9)), (K9,K(9)), (L9,L(9))
C     EQUIVALENCE (I10,I(10)), (K10,K(10)), (L10,L(10))
C     EQUIVALENCE (K11,K(11))
C     .. Scalar Arguments ..
      INTEGER           PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
      INTEGER           SYM(21)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           I1, I10, I2, I3, I4, I5, I6, I7, I8, I9, J, JJ,
     *                  K1, K10, K11, K2, K3, K4, K5, K6, K7, K8, K9,
     *                  KK, L1, L10, L2, L3, L4, L5, L6, L7, L8, L9,
     *                  LEVEL, LOOP, NEST
C     .. Local Arrays ..
      INTEGER           I(20), K(20), L(20)
C     .. Data statements ..
      DATA              NEST/20/
      DATA              LOOP/10/
C     .. Executable Statements ..
      IF (SYM(1).EQ.0) GO TO 360
      DO 20 J = 1, NEST
         L(J) = 1
         I(J) = 1
   20 CONTINUE
      KK = PTS
      DO 40 J = 1, NEST
         IF (SYM(J).EQ.0) GO TO 60
         L(J) = KK
         I(J) = KK/SYM(J)
         KK = KK/SYM(J)
   40 CONTINUE
   60 CONTINUE
C
      L1 = L(1)
      L2 = L(2)
      L3 = L(3)
      L4 = L(4)
      L5 = L(5)
      L6 = L(6)
      L7 = L(7)
      L8 = L(8)
      L9 = L(9)
      L10 = L(10)
      I1 = I(1)
      I2 = I(2)
      I3 = I(3)
      I4 = I(4)
      I5 = I(5)
      I6 = I(6)
      I7 = I(7)
      I8 = I(8)
      I9 = I(9)
      I10 = I(10)
C
      KK = 0
      LEVEL = NEST
      K(LEVEL) = 1
      GO TO 100
   80 CONTINUE
      IF (LEVEL.GE.NEST) GO TO 360
      LEVEL = LEVEL + 1
      K(LEVEL) = K(LEVEL) + I(LEVEL)
      IF (K(LEVEL).GT.L(LEVEL)) GO TO 80
  100 CONTINUE
      LEVEL = LEVEL - 1
      DO 120 J = LOOP, LEVEL
         JJ = LEVEL + LOOP - J
         K(JJ) = K(JJ+1)
  120 CONTINUE
      K11 = K(11)
      DO 340 K10 = K11, L10, I10
         DO 320 K9 = K10, L9, I9
            DO 300 K8 = K9, L8, I8
               DO 280 K7 = K8, L7, I7
                  DO 260 K6 = K7, L6, I6
                     DO 240 K5 = K6, L5, I5
                        DO 220 K4 = K5, L4, I4
                           DO 200 K3 = K4, L3, I3
                              DO 180 K2 = K3, L2, I2
                                 DO 160 K1 = K2, L1, I1
                                    KK = KK + 1
                                    IF (KK.GE.K1) GO TO 140
                                    T = X(KK)
                                    X(KK) = X(K1)
                                    X(K1) = T
                                    T = Y(KK)
                                    Y(KK) = Y(K1)
                                    Y(K1) = T
  140                               CONTINUE
  160                            CONTINUE
  180                         CONTINUE
  200                      CONTINUE
  220                   CONTINUE
  240                CONTINUE
  260             CONTINUE
  280          CONTINUE
  300       CONTINUE
  320    CONTINUE
  340 CONTINUE
      LEVEL = LOOP
      GO TO 80
C
  360 CONTINUE
      RETURN
      END
